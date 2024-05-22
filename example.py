import jax
import numpy as np
import netket as nk
from netket.operator.spin import sigmax, sigmaz
#

N = 10
hi = nk.hilbert.Spin(s=0.5, N=10)
Gamma = -1.0
H = sum([Gamma * sigmax(hi, i) * sigmax(hi, (i + 1) % N) for i in range(N)])
H += sum([Gamma * sigmax(hi, i) for i in range(N)])

# H.apply(hi.random_state(jax.random.PRNGKey(1234), N))

x = np.ones(hi.size)
# x[0] = 1.0

print(H.get_conn(x))
(conn, elems) = H.get_conn(x)

print(conn)
print(conn.shape)
# print(hi.size)
