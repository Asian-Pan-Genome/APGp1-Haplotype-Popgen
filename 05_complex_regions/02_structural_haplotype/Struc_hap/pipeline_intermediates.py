#!/usr/bin/env python3
"""Bounded cleanup of disposable pipeline intermediates.

Callers register only files created by the current successful operation. Merely
having a familiar filename is never evidence of ownership: these helpers do not
scan directories, recurse, remove directories, or delete unregistered files.
Reports, metadata, mapping evidence, snapshots, and logs are not disposable.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import stat
import uuid
from pathlib import Path


_MAFFT_FASTA = re.compile(r"(?:source|compare)_[0-9]+(?:\.aln)?\.fa\Z")
_SELF_FASTA = re.compile(r"iterations/iteration_[0-9]{4,}/bundle(?:_after)?\.fa\Z")
_SELF_DIRECTORY = re.compile(r"iteration_[0-9]{4,}\Z")
_RANK_FASTAS = {"rank_population.fa", "remap/split_evidence/previous_ranks.fa"}
_FINAL_DUPLICATE = "mapping_split/bundle.fa"
_FINAL_COPY = "final_bundle.fa"
_FINAL_SUMMARY = "mapping_split/mapping_split_summary.json"


def _empty_report(keep_intermediates):
    return {"deleted": [], "missing": [], "skipped": [],
            "deleted_count": 0, "deleted_bytes": 0,
            "keep_intermediates": bool(keep_intermediates)}


def _registered_relative_path(stage, path):
    """Normalize a registered exact path without following filesystem links."""
    value = Path(path)
    if ".." in value.parts:
        raise ValueError("parent_traversal")
    boundaries = (stage, Path(os.path.realpath(stage.parent)) / stage.name)
    if value.is_absolute():
        for boundary in boundaries:
            try:
                value = value.relative_to(boundary)
                break
            except ValueError:
                continue
        else:
            raise ValueError("outside_stage")
    else:
        # Pipeline code also naturally passes directory / generated_filename
        # when its output directory is relative to the working directory.
        for boundary in boundaries:
            try:
                value = Path(os.path.abspath(value)).relative_to(boundary)
                break
            except ValueError:
                continue
    if not value.parts or value == Path("."):
        raise ValueError("not_a_file_path")
    return value.as_posix()


def _whitelisted(stage, relative):
    if "/" not in relative and _MAFFT_FASTA.fullmatch(relative):
        return True
    if relative in _RANK_FASTAS or relative == _FINAL_DUPLICATE:
        return True
    if _SELF_FASTA.fullmatch(relative):
        return True
    return (stage.parent.name == "iterations" and _SELF_DIRECTORY.fullmatch(stage.name)
            and relative in {"bundle.fa", "bundle_after.fa"})


def _directory_flags():
    return os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW


def _open_stage(stage):
    # Ancestors outside the selected boundary may include platform aliases such
    # as macOS /var -> /private/var. The boundary itself and every parent below
    # it must be real directories; descriptor-relative access prevents escape.
    info = stage.lstat()
    if not stat.S_ISDIR(info.st_mode):
        raise ValueError("stage_not_regular_directory")
    descriptor = os.open(stage, _directory_flags())
    opened = os.fstat(descriptor)
    if (opened.st_dev, opened.st_ino) != (info.st_dev, info.st_ino):
        os.close(descriptor)
        raise ValueError("stage_changed")
    return descriptor


def _open_parent(stage_descriptor, relative):
    """Open only exact registered parents, refusing all links below the stage."""
    parts = Path(relative).parts
    descriptor = os.dup(stage_descriptor)
    try:
        for part in parts[:-1]:
            replacement = os.open(part, _directory_flags(), dir_fd=descriptor)
            os.close(descriptor)
            descriptor = replacement
        return descriptor, parts[-1]
    except BaseException:
        os.close(descriptor)
        raise


def _protected_paths(stage, paths):
    lexical, identities = set(), set()
    for path in paths:
        value = Path(path)
        # Protected paths are input paths, relative to the caller's working
        # directory, as are the command-line paths elsewhere in the pipeline.
        absolute = Path(os.path.abspath(value))
        lexical.add(absolute)
        if not value.is_absolute():
            lexical.add(stage / value)
        try:
            info = value.stat()
        except OSError:
            continue
        identities.add((info.st_dev, info.st_ino))
    return lexical, identities


def _same_bytes(stage_descriptor, relative, candidate_info):
    """Verify the final published FASTA is an exact, safe redundant copy."""
    descriptors = []
    try:
        for path in (relative, _FINAL_COPY):
            parent, name = _open_parent(stage_descriptor, path)
            try:
                info = os.stat(name, dir_fd=parent, follow_symlinks=False)
                if not stat.S_ISREG(info.st_mode):
                    return False
                descriptor = os.open(name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=parent)
                descriptors.append(descriptor)
                opened = os.fstat(descriptor)
                if ((opened.st_dev, opened.st_ino, opened.st_size)
                        != (info.st_dev, info.st_ino, info.st_size)):
                    return False
            finally:
                os.close(parent)
        left, right = (os.fstat(descriptor) for descriptor in descriptors)
        if ((left.st_dev, left.st_ino, left.st_size)
                != (candidate_info.st_dev, candidate_info.st_ino, candidate_info.st_size)
                or left.st_size != right.st_size):
            return False
        while True:
            first, second = (os.read(descriptor, 1024 * 1024) for descriptor in descriptors)
            if first != second:
                return False
            if not first:
                return True
    except (OSError, ValueError):
        return False
    finally:
        for descriptor in descriptors:
            os.close(descriptor)


def _alias_final_summary(stage_descriptor, stage, registered, protected, identities):
    """Atomically retain a valid summary reference before removing its duplicate.

    No existing metadata is edited without independent current-run registration.
    If metadata cannot safely be rewritten, the caller retains the FASTA copy.
    """
    parent = descriptor = temporary_descriptor = None
    temporary_name = None
    try:
        parent, name = _open_parent(stage_descriptor, _FINAL_SUMMARY)
        try:
            info = os.stat(name, dir_fd=parent, follow_symlinks=False)
        except FileNotFoundError:
            return "not_present"
        if (not registered or not stat.S_ISREG(info.st_mode)
                or stage / _FINAL_SUMMARY in protected
                or (info.st_dev, info.st_ino) in identities):
            return None
        descriptor = os.open(name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=parent)
        opened = os.fstat(descriptor)
        if (opened.st_dev, opened.st_ino) != (info.st_dev, info.st_ino):
            return None
        pieces = []
        while True:
            piece = os.read(descriptor, 1024 * 1024)
            if not piece:
                break
            pieces.append(piece)
        summary = json.loads(b"".join(pieces))
        outputs = summary.get("outputs") if isinstance(summary, dict) else None
        if not isinstance(outputs, dict) or not isinstance(outputs.get("bundle.fa"), str):
            return None
        reference = _registered_relative_path(stage, outputs["bundle.fa"])
        if reference == _FINAL_COPY:
            return "already_aliased"
        if reference != _FINAL_DUPLICATE:
            return None
        retained = str(Path(os.path.realpath(stage.parent)) / stage.name / _FINAL_COPY)
        outputs["bundle.fa"] = retained
        summary["intermediate_retention"] = {
            "bundle_fasta_alias": {"redundant_copy": _FINAL_DUPLICATE,
                                   "retained_copy": retained},
            "policy": "registered redundant FASTA copy; reports and evidence retained",
        }
        content = (json.dumps(summary, indent=2, sort_keys=True) + "\n").encode("utf-8")
        temporary_name = f".pipeline_intermediates-summary-{uuid.uuid4().hex}.tmp"
        temporary_descriptor = os.open(temporary_name,
                                       os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                                       stat.S_IMODE(info.st_mode), dir_fd=parent)
        offset = 0
        while offset < len(content):
            offset += os.write(temporary_descriptor, content[offset:])
        os.fsync(temporary_descriptor)
        current = os.stat(name, dir_fd=parent, follow_symlinks=False)
        if ((current.st_dev, current.st_ino, current.st_size, current.st_mtime_ns)
                != (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns)
                or not stat.S_ISREG(current.st_mode)):
            return None
        os.replace(temporary_name, name, src_dir_fd=parent, dst_dir_fd=parent)
        temporary_name = None
        return "aliased"
    except FileNotFoundError:
        # No summary exists, hence no summary output can become dangling.
        return "not_present"
    except (OSError, ValueError, TypeError, UnicodeError):
        return None
    finally:
        for value in (descriptor, temporary_descriptor):
            if value is not None:
                os.close(value)
        if temporary_name is not None and parent is not None:
            try:
                os.unlink(temporary_name, dir_fd=parent)
            except OSError:
                pass
        if parent is not None:
            os.close(parent)


def cleanup_generated_files(directory, generated_paths, protected_paths=(),
                            keep_intermediates=False):
    """Unlink registered, whitelisted regular files within one stage boundary.

    ``generated_paths`` must contain exact paths newly written by the current
    operation, not filenames discovered in an existing directory. Both absolute
    paths under ``directory`` and paths relative to it are accepted. Call after
    successful consumption/validation, and protect original inputs plus retained
    first/last self-alignment FASTAs. A failed operation must not call cleanup.
    Symlinks, unknown names, protected files (including hardlinks), and differing
    final-refinement copies are retained and recorded with reasons.
    """
    stage = Path(os.path.abspath(directory))
    report = _empty_report(keep_intermediates)
    protected, identities = _protected_paths(stage, protected_paths)
    candidates, seen, generated_summary = [], set(), False
    for path in generated_paths:
        try:
            relative = _registered_relative_path(stage, path)
        except ValueError as error:
            report["skipped"].append({"path": str(path), "reason": str(error)})
            continue
        if relative in seen:
            continue
        seen.add(relative)
        if relative == _FINAL_SUMMARY:
            generated_summary = True
            report["skipped"].append({"path": relative, "reason": "retained_metadata"})
        elif not _whitelisted(stage, relative):
            report["skipped"].append({"path": relative, "reason": "not_whitelisted"})
        elif stage / relative in protected:
            report["skipped"].append({"path": relative, "reason": "protected_path"})
        elif keep_intermediates:
            report["skipped"].append({"path": relative, "reason": "keep_intermediates"})
        else:
            candidates.append(relative)
    if not candidates:
        return report
    try:
        stage_descriptor = _open_stage(stage)
    except FileNotFoundError:
        report["missing"].extend(candidates)
        return report
    except (OSError, ValueError):
        report["skipped"].extend({"path": path, "reason": "unsafe_stage_directory"}
                                 for path in candidates)
        return report
    try:
        for relative in candidates:
            parent = None
            try:
                parent, name = _open_parent(stage_descriptor, relative)
                info = os.stat(name, dir_fd=parent, follow_symlinks=False)
                if not stat.S_ISREG(info.st_mode):
                    reason = "not_regular_file"
                elif (info.st_dev, info.st_ino) in identities:
                    reason = "protected_inode"
                elif relative == _FINAL_DUPLICATE and not _same_bytes(stage_descriptor, relative, info):
                    reason = "final_copy_not_identical"
                else:
                    if relative == _FINAL_DUPLICATE:
                        alias = _alias_final_summary(stage_descriptor, stage, generated_summary,
                                                     protected, identities)
                        if alias is None:
                            report["skipped"].append({"path": relative,
                                                       "reason": "metadata_alias_unavailable"})
                            continue
                        if alias in {"aliased", "already_aliased"}:
                            report.setdefault("aliased_outputs", []).append({
                                "metadata": _FINAL_SUMMARY, "output": "bundle.fa",
                                "retained_path": str(Path(os.path.realpath(stage.parent))
                                                     / stage.name / _FINAL_COPY)})
                    # Recheck identity immediately before unlink. No recursive
                    # or path-following removal is ever used.
                    current = os.stat(name, dir_fd=parent, follow_symlinks=False)
                    if ((current.st_dev, current.st_ino, current.st_size, current.st_mtime_ns)
                            != (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns)
                            or not stat.S_ISREG(current.st_mode)):
                        reason = "file_changed"
                    else:
                        os.unlink(name, dir_fd=parent)
                        report["deleted"].append({"path": relative, "bytes": info.st_size})
                        report["deleted_count"] += 1
                        report["deleted_bytes"] += info.st_size
                        continue
                report["skipped"].append({"path": relative, "reason": reason})
            except FileNotFoundError:
                report["missing"].append(relative)
            except OSError:
                report["skipped"].append({"path": relative, "reason": "unsafe_or_unremovable_path"})
            finally:
                if parent is not None:
                    os.close(parent)
    finally:
        os.close(stage_descriptor)
    return report


def cleanup_rank_round(directory, protected_paths=(), keep_intermediates=False, *,
                       generated_paths=()):
    """Reduce one pipeline-owned rank directory to raw PAFs and additions.

    Only exact top-level files created by the rank driver and files inside its
    three pipeline-owned subdirectories are eligible. Unknown top-level files
    are deliberately ignored. Raw ``*.paf`` files, ``additions.tsv``, the round
    manifest, and the cleanup audit are retained. Symlinks and protected input
    identities are never unlinked.
    """
    stage = Path(os.path.abspath(directory))
    report = (cleanup_generated_files(directory, generated_paths, protected_paths,
                                      keep_intermediates)
              if generated_paths else _empty_report(keep_intermediates))
    protected, identities = _protected_paths(stage, protected_paths)
    retained = {"additions.tsv", "round_artifacts.json", "intermediate_cleanup.json"}
    top_level = {
        "old_library_coverage.bed", "skipped_candidates.tsv", "source_bundles.pgr.bed",
        "library_initial.fa", "library_initial.metadata.bed",
        "library_before.fa", "library_before.metadata.bed",
        "library_before_self.fa", "library_before_self.metadata.bed",
        "library_after.fa", "library_after.metadata.bed", "rank_population.fa",
    }
    candidates = []
    try:
        children = list(stage.iterdir())
    except FileNotFoundError:
        return report
    for child in children:
        if child.name in {"trf", "remap", "self_alignment"}:
            if child.is_symlink() or not child.is_dir():
                report["skipped"].append({"path": child.name, "reason": "unsafe_owned_directory"})
                continue
            for path in child.rglob("*"):
                if path.is_symlink():
                    report["skipped"].append({"path": str(path.relative_to(stage)), "reason": "symlink"})
                elif path.is_file() and path.suffix != ".paf":
                    candidates.append(path)
        elif (child.name in top_level
              or re.fullmatch(r"(?:source|compare)_[0-9]+\.mafft\.log", child.name)) \
                and child.name not in retained:
            candidates.append(child)
    for path in candidates:
        relative = str(path.relative_to(stage))
        if path in protected:
            report["skipped"].append({"path": relative, "reason": "protected_path"})
            continue
        try:
            info = path.lstat()
        except FileNotFoundError:
            report["missing"].append({"path": relative})
            continue
        if not stat.S_ISREG(info.st_mode):
            report["skipped"].append({"path": relative, "reason": "not_regular_file"})
            continue
        if (info.st_dev, info.st_ino) in identities:
            report["skipped"].append({"path": relative, "reason": "protected_identity"})
            continue
        if keep_intermediates:
            report["skipped"].append({"path": relative, "reason": "keep_intermediates"})
            continue
        path.unlink()
        report["deleted"].append({"path": relative, "bytes": info.st_size})
        report["deleted_count"] += 1
        report["deleted_bytes"] += info.st_size
    return report


def cleanup_self_iteration(directory, iteration, final_iteration=False,
                           protected_paths=(), keep_intermediates=False, *,
                           generated_paths=()):
    """Prune a previous iteration after the current input was safely written.

    The first input and current (possibly final) input/output remain protected.
    Only the previous iteration's caller-registered FASTAs can be removed.
    ``final_iteration`` documents the lifecycle boundary; it never broadens it.
    """
    if isinstance(iteration, bool) or not isinstance(iteration, int) or iteration < 1:
        raise ValueError("iteration must be a positive integer")
    stage = Path(os.path.abspath(directory))
    previous = f"iterations/iteration_{iteration - 1:04d}"
    allowed = {f"{previous}/bundle_after.fa"}
    if iteration > 2:
        allowed.add(f"{previous}/bundle.fa")
    registered, retained = [], []
    for path in generated_paths:
        try:
            relative = _registered_relative_path(stage, path)
        except ValueError as error:
            retained.append({"path": str(path), "reason": str(error)})
            continue
        if iteration > 1 and relative in allowed:
            registered.append(path)
        else:
            retained.append({"path": relative, "reason": "retained_self_boundary"})
    current = stage / "iterations" / f"iteration_{iteration:04d}"
    protected = [*protected_paths, stage / "iterations/iteration_0001/bundle.fa",
                 current / "bundle.fa", current / "bundle_after.fa"]
    report = cleanup_generated_files(stage, registered, protected, keep_intermediates)
    report["skipped"].extend(retained)
    return report


def cleanup_final_refinement(directory, protected_paths=(), keep_intermediates=False, *,
                             generated_paths=()):
    """Remove only a registered exact redundant mapping-split FASTA copy."""
    stage = Path(os.path.abspath(directory))
    registered, retained = [], []
    for path in generated_paths:
        try:
            relative = _registered_relative_path(stage, path)
            if relative in {_FINAL_DUPLICATE, _FINAL_SUMMARY}:
                registered.append(path)
            else:
                retained.append({"path": relative, "reason": "not_disposable_for_operation"})
        except ValueError:
            # Delegate invalid registrations to the shared safety audit.
            registered.append(path)
    report = cleanup_generated_files(stage, registered, protected_paths, keep_intermediates)
    report["skipped"].extend(retained)
    return report


def _write_audit(directory, path, output):
    """Write the explicit audit output without following directory/file links."""
    stage = Path(os.path.abspath(directory))
    relative = _registered_relative_path(stage, path)
    if relative != "intermediate_cleanup.json":
        raise ValueError("audit path must be intermediate_cleanup.json in the stage directory")
    stage_descriptor = _open_stage(stage)
    temporary_descriptor = None
    temporary_name = None
    try:
        try:
            info = os.stat(relative, dir_fd=stage_descriptor, follow_symlinks=False)
            if not stat.S_ISREG(info.st_mode):
                raise ValueError("audit output is not a regular file")
        except FileNotFoundError:
            pass
        temporary_name = f".pipeline_intermediates-audit-{uuid.uuid4().hex}.tmp"
        temporary_descriptor = os.open(temporary_name,
                                       os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                                       0o644, dir_fd=stage_descriptor)
        content, offset = output.encode("utf-8"), 0
        while offset < len(content):
            offset += os.write(temporary_descriptor, content[offset:])
        os.fsync(temporary_descriptor)
        os.replace(temporary_name, relative,
                   src_dir_fd=stage_descriptor, dst_dir_fd=stage_descriptor)
        temporary_name = None
    finally:
        if temporary_descriptor is not None:
            os.close(temporary_descriptor)
        if temporary_name is not None:
            try:
                os.unlink(temporary_name, dir_fd=stage_descriptor)
            except OSError:
                pass
        os.close(stage_descriptor)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=("final-refinement",))
    parser.add_argument("directory", type=Path)
    parser.add_argument("--generated-path", type=Path, action="append", default=[])
    parser.add_argument("--protected-path", type=Path, action="append", default=[])
    parser.add_argument("--keep-intermediates", action="store_true")
    parser.add_argument("--audit-path", type=Path)
    args = parser.parse_args(argv)
    report = cleanup_final_refinement(args.directory, args.protected_path,
                                      args.keep_intermediates, generated_paths=args.generated_path)
    report["operation"] = args.operation
    output = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.audit_path is not None:
        # This is an explicitly requested audit output, never a cleanup target.
        _write_audit(args.directory, args.audit_path, output)
    else:
        print(output, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
