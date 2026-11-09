.PHONY: test lint check-data reproduce

test:
	uv run python -m unittest discover -s tests -v

lint:
	uv run ruff check src tests
	uv run ruff format --check src tests

check-data:
	cd data && shasum -a 256 -c SHA256SUMS

reproduce: check-data
	uv run gene-network-analysis run
