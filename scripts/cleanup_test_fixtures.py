#!/usr/bin/env python3
"""Compatibility entry point for the guarded test-fixture cleanup service.

The cleanup implementation lives in ``backend.cleanup_test_fixtures`` so ORM
and raw SQLite verification share the explicitly configured TEST/E2E database
instead of reopening the repository production database by relative path.
"""

from backend.cleanup_test_fixtures import run_cleanup


def main() -> dict:
    return run_cleanup()


if __name__ == "__main__":
    main()
