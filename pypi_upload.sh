rm -r build
rm -r dist
python setup.py sdist
python setup.py bdist_wheel
twine upload dist/*