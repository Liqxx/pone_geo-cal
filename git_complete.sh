#!/bin/bash

note="$1"

git add -A
git commit -m "$note"
git push