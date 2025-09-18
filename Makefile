.PHONY: help install install-dev test lint format clean build upload
.DEFAULT_GOAL := help

help: ## Show this help message
	@echo 'Usage: make [target]'
	@echo ''
	@echo 'Targets:'
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "  %-15s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

install: ## Install the package
	pip install -e .

install-dev: ## Install the package in development mode with dev dependencies
	pip install -e ".[dev]"

test: ## Run tests
	pytest

test-cov: ## Run tests with coverage report
	pytest --cov=detaxa --cov-report=html --cov-report=term

lint: ## Run code linting
	flake8 src/detaxa
	mypy src/detaxa

format: ## Format code with black
	black src/detaxa

format-check: ## Check code formatting without making changes
	black --check src/detaxa

clean: ## Clean build artifacts
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .pytest_cache/
	rm -rf htmlcov/
	find . -type d -name __pycache__ -delete
	find . -type f -name "*.pyc" -delete

build: clean ## Build the package
	python -m build

upload-test: build ## Upload to TestPyPI
	python -m twine upload --repository testpypi dist/*

upload: build ## Upload to PyPI
	python -m twine upload dist/*

check: ## Check package metadata
	python -m twine check dist/*
