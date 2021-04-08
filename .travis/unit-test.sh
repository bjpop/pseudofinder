#!/bin/bash

set -e
errors=0

# Run unit tests
python psuedofinder/psuedofinder_test.py || {
    echo "'python python/psuedofinder/psuedofinder_test.py' failed"
    let errors+=1
}

# Check program style
pylint -E psuedofinder/*.py || {
    echo 'pylint -E psuedofinder/*.py failed'
    let errors+=1
}

[ "$errors" -gt 0 ] && {
    echo "There were $errors errors found"
    exit 1
}

echo "Ok : Python specific tests"
