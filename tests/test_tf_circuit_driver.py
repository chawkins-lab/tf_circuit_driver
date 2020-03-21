#!/usr/bin/env python

"""Tests for `tf_circuit_driver` package."""


import unittest
from click.testing import CliRunner

from tf_circuit_driver import tf_circuit_driver
from tf_circuit_driver import cli


class TestTf_circuit_driver(unittest.TestCase):
    """Tests for `tf_circuit_driver` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_000_something(self):
        """Test something."""

    def test_command_line_interface(self):
        """Test the CLI."""
        runner = CliRunner()
        result = runner.invoke(cli.main)
        assert result.exit_code == 0
        assert 'tf_circuit_driver.cli.main' in result.output
        help_result = runner.invoke(cli.main, ['--help'])
        assert help_result.exit_code == 0
        assert '--help  Show this message and exit.' in help_result.output
