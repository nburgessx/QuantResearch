python pacman.py
python pacman.py -p ReflexAgent
python pacman.py -p ReflexAgent -l testClassic
python pacman.py -p ReflexAgent -l testClassic
python pacman.py --frameTime 0 -p ReflexAgent -k 1
python pacman.py --frameTime 0 -p ReflexAgent -k 2
python pacman.py -p ReflexAgent -l openClassic -n 10 -q
python pacman.py -p MinimaxAgent -l minimaxClassic -a depth=4
python pacman.py -p MinimaxAgent -l trappedClassic -a depth=3
python pacman.py -p AlphaBetaAgent -a depth=3 -l smallClassic
python pacman.py -p AlphaBetaAgent -l trappedClassic -a depth=3 -q -n 10
python pacman.py -p ExpectimaxAgent -l trappedClassic -a depth=3 -q -n 10
python pacman.py -l smallClassic -p ExpectimaxAgent -a evalFn=better -q -n 10
python pacman.py -l contestClassic -p ContestAgent -g DirectionalGhost -q -n 10
