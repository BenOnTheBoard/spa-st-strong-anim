from bruteforce import STSMBruteForce
from SPAST_Strong_Anim import SPAST_STRONG_ANIM

filename = "examples/problematic/K55.txt"
bruteforcer = STSMBruteForce(filename)
result = bruteforcer.choose(1)
print(result)

strong_solver = SPAST_STRONG_ANIM(filename)
strong_solver.inner_repeat()
print(strong_solver.G)
