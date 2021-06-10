#!/usr/bin/env python
# coding: utf-8

import json
import typer

from pathlib import Path


app = typer.Typer()


@app.command()
def main(nextflow_schema_json: Path):
    """Convert nextflow_schema.json to Markdown.

    Outputs Markdown to stdout.

    Usage:

    $ ./nextflow_schema_json_to_markdown.py nextflow_schema.json > workflow_params.md
    """
    with open(nextflow_schema_json) as fh:
        nextflow_schema = json.load(fh)

    defs = nextflow_schema['definitions']

    mdl = []
    for def_key, def_dict in defs.items():
        def_title = def_dict.get('title', None)
        if def_title is None:
            raise ValueError(f'Nextflow schema definition "{def_key}" does not have a `title`! Please check your Nextflow schema JSON!')
        mdl.append(f'### {def_title}\n')
        def_desc = def_dict.get('description', None)
        if def_desc:
            mdl.append(f'{def_desc}\n')
        props = def_dict.get('properties', {})
        required = set(def_dict.get('required', []))
        for param, param_info in props.items():
            param_type = param_info.get('type', None)
            param_default = param_info.get('default', None)
            param_min = param_info.get('minimum', None)
            param_max = param_info.get('maximum', None)
            param_desc = param_info.get('description', None)
            param_help_text = param_info.get('help_text', None)
            is_req = param in required
            mdl.append(f'#### `--{param}`\n')
            if is_req:
                mdl.append(f'- **Required**')
            else:
                mdl.append(f'- Optional')
            if param_type:
                mdl.append(f'- Type: {param_type}')
            if param_min and param_max:
                mdl.append(f'  - Range: {param_min} - {param_max}')
            elif param_min:
                mdl.append(f'  - Minimum: {param_min}')
            elif param_max:
                mdl.append(f'  - Maximum: {param_max}')
            if param_default:
                mdl.append(f'- Default: `{param_default}`')
            mdl.append(f'\n{param_desc}')
            if param_help_text:
                mdl.append(f'\n> **NOTE:** {param_help_text}')
            mdl.append('\n')
    print('\n'.join(mdl))


if __name__ == '__main__':
    app()
