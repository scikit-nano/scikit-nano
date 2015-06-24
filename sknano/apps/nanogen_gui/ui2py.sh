#!/usr/bin/env bash

pyuic5 nanogen.ui > _ui_nanogen.py
pyuic5 nanogen_swnts.ui > _ui_nanogen_swnts.py
pyuic5 nanogen_mwnts.ui > _ui_nanogen_mwnts.py
pyuic5 nanogen_graphene.ui > _ui_nanogen_graphene.py
pyuic5 nanogen_fullerenes.ui > _ui_nanogen_fullerenes.py
pyuic5 nanogen_bulk_structures.ui > _ui_nanogen_bulk_structures.py
