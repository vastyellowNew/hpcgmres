if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    sudo apt-get update -y
    sudo apt-get install -y lcov
    sudo apt-get install -y valgrind
else if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    brew update
    brew install lcov
    fi
fi