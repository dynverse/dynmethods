bucket="travis-python-archives"
vers="3.6"
lang="python"
ext="bz2"
PYENV_PATH_FILE="/etc/profile.d/pyenv.sh"
archive_basename="${lang}-${vers}"
archive_filename="${archive_basename}.tar.bz2"
travis_host_os=$(lsb_release -is | tr 'A-Z' 'a-z')
travis_rel_version=$(lsb_release -rs)
archive_url=https://s3.amazonaws.com/${bucket}/binaries/${travis_host_os}/${travis_rel_version}/$(uname -m)/${archive_filename}

echo "Downloading archive: ${archive_url}"
curl -sSf -o ${archive_filename} ${archive_url}
sudo tar xjf ${archive_filename} --directory /
rm ${archive_filename}

echo 'export PATH=/opt/python/${vers}/bin:$PATH' | sudo tee -a ${PYENV_PATH_FILE} &>/dev/null
export PATH=/opt/python/${vers}/bin:$PATH
