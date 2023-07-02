

show:
	@echo "Current environment:"
	@poetry env info

install:
	@poetry install

format-check:
	@black --check --diff .

format:
	@black .
