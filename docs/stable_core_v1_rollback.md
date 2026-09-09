# Stable Core v1 database rollback

The authoritative pre-migration backup is recorded in
`validation/stable_core_v1_pre_migration_manifest.json`. It is intentionally
stored outside Git.

Rollback is a controlled disaster-recovery operation, not an application
migration:

1. stop the existing `drug-opt.service`;
2. verify the backup SHA-256 against the pre-migration manifest;
3. make a separate backup of the post-migration database;
4. replace `drug_opt.db` with the verified pre-migration copy;
5. start the existing service and verify physical and business integrity.

Because restoring the pre-migration file discards transactions made after the
migration, rollback requires an explicit recovery decision. Git rollback alone
does not roll back SQLite state.
