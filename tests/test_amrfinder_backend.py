import pytest
import subprocess
from unittest.mock import patch
from genome_firewall.run_amrfinder import select_backend

@patch('platform.system')
@patch('shutil.which')
@patch('subprocess.run')
def test_select_backend_auto_windows(mock_run, mock_which, mock_system):
    mock_system.return_value = 'Windows'
    # Windows always selects docker for auto
    assert select_backend('auto') == 'docker'
    
@patch('platform.system')
@patch('shutil.which')
def test_select_backend_auto_linux_with_native(mock_which, mock_system):
    mock_system.return_value = 'Linux'
    mock_which.return_value = '/usr/bin/amrfinder'
    assert select_backend('auto') == 'native'

@patch('platform.system')
@patch('shutil.which')
def test_select_backend_auto_linux_without_native(mock_which, mock_system):
    mock_system.return_value = 'Linux'
    mock_which.return_value = None  # no amrfinder
    assert select_backend('auto') == 'docker'

@patch('shutil.which')
def test_select_backend_native_missing(mock_which):
    mock_which.return_value = None
    with pytest.raises(RuntimeError, match="not in PATH"):
        select_backend('native')

@patch('shutil.which')
@patch('subprocess.run')
def test_select_backend_docker_running(mock_run, mock_which):
    mock_which.return_value = '/usr/bin/docker'
    mock_run.return_value.returncode = 0
    assert select_backend('docker') == 'docker'

@patch('shutil.which')
@patch('subprocess.run')
def test_select_backend_docker_not_running(mock_run, mock_which):
    mock_which.return_value = '/usr/bin/docker'
    mock_run.side_effect = subprocess.CalledProcessError(1, 'docker info')
    with pytest.raises(RuntimeError, match="daemon is not running"):
        select_backend('docker')
