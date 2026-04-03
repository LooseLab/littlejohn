# `robin password`

Configure the **NiceGUI** workflow monitor login password.

For the full startup sequence (disclaimer, when the GUI prompts for a password), see **[What happens at startup](../getting-started/startup.md)**.

## `robin password set`

```bash
robin password set
```

- Prompts twice for the new password (no echo).
- Replaces any previously stored password.

If the GUI password module fails to import (minimal install / broken env), the command prints an error and exits non-zero.

## When it applies

The password is used when the **[`robin workflow`](workflow.md)** GUI is enabled (`--with-gui`, default) and password protection is active in the GUI layer — see application logs and on-screen prompts if login is required.

## Related

- [`robin workflow` GUI options](workflow.md#gui-nicegui)  
- [CLI overview](index.md)  
