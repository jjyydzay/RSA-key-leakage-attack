from __future__ import print_function
import logging
import time
import os
from fpylll import IntegerMatrix, LLL


debug = True
strict = False
helpful_only = True
dimension_min = 7 # stop removing if lattice reaches that dimension

def reduce_lattice(B, delta=0.8):
    verbose = False
    #Reduce basis
    start = time.time() #Can't use cputime here, because it does not catch os.system(...).  
    if verbose:
        print("Reducing a ",B.nrows()," x ",B.ncols()," lattice...")

    with open("basis.tmp", "w+") as f:
    #Would prefer to simply use IntegerMatrix here, but IntegerMatrix.__str()__ results in crash.
        B_str = B.str()
        B_str = '\n'.join(' '.join(line.split()) for line in B_str.split('\n'))
        f.write( "[\n" + B_str + "\n]")
    # print("before success")
    success = os.system("flatter -v basis.tmp basis_out.tmp >/dev/null 2>&1")
    # success = os.system("flatter -v basis.tmp basis_out.tmp")

    # print(success)
    os.remove("basis.tmp")

    if success == 0:
        B_LLL = matrix(IntegerMatrix.from_file("basis_out.tmp"))
        os.remove("basis_out.tmp")
    else:
        if verbose:
            print("flatter not found. Resorting to FPLLL.")
        B_LLL = B.LLL(delta)

    stop = time.time()
    if verbose:
        print("Finished basis Flatter reduction. Time: %fs." % (stop-start))
    return B_LLL

# display stats on helpful vectors
def helpful_vectors(BB, modulus):
    nothelpful = 0
    for ii in range(BB.dimensions()[0]):
        if BB[ii,ii] >= modulus:
            nothelpful += 1
 
    print(nothelpful, "/", BB.dimensions()[0], " vectors are not helpful")
 
# display matrix picture with 0 and X
def matrix_overview(BB, bound):
    for ii in range(BB.dimensions()[0]):
        a = ('%02d ' % ii)
        for jj in range(BB.dimensions()[1]):
            a += '0' if BB[ii,jj] == 0 else 'X'
            if BB.dimensions()[0] < 60:
                a += ' '
        if BB[ii, ii] >= bound:
            a += '~'
        print(a)
 
# tries to remove unhelpful vectors
# we start at current = n-1 (last vector)
def remove_unhelpful(BB, monomials, bound, current):
    # end of our recursive function
    if current == -1 or BB.dimensions()[0] <= dimension_min:
        return BB
 
    # we start by checking from the end
    for ii in range(current, -1, -1):
        # if it is unhelpful:
        if BB[ii, ii] >= bound:
            affected_vectors = 0
            affected_vector_index = 0
            # let's check if it affects other vectors
            for jj in range(ii + 1, BB.dimensions()[0]):
                # if another vector is affected:
                # we increase the count
                if BB[jj, ii] != 0:
                    affected_vectors += 1
                    affected_vector_index = jj
 
            # level:0
            # if no other vectors end up affected
            # we remove it
            if affected_vectors == 0:
                print("* removing unhelpful vector", ii)
                BB = BB.delete_columns([ii])
                BB = BB.delete_rows([ii])
                monomials.pop(ii)
                BB = remove_unhelpful(BB, monomials, bound, ii-1)
                return BB
 
            # level:1
            # if just one was affected we check
            # if it is affecting someone else
            elif affected_vectors == 1:
                affected_deeper = True
                for kk in range(affected_vector_index + 1, BB.dimensions()[0]):
                    # if it is affecting even one vector
                    # we give up on this one
                    if BB[kk, affected_vector_index] != 0:
                        affected_deeper = False
                # remove both it if no other vector was affected and
                # this helpful vector is not helpful enough
                # compared to our unhelpful one
                if affected_deeper and abs(bound - BB[affected_vector_index, affected_vector_index]) < abs(bound - BB[ii, ii]):
                    print("* removing unhelpful vectors", ii, "and", affected_vector_index)
                    BB = BB.delete_columns([affected_vector_index, ii])
                    BB = BB.delete_rows([affected_vector_index, ii])
                    monomials.pop(affected_vector_index)
                    monomials.pop(ii)
                    BB = remove_unhelpful(BB, monomials, bound, ii-1)
                    return BB
    # nothing happened
    return BB
 
def boneh_durfee(pol, modulus, mm, tt, XX, YY):
 
    # substitution (Herrman and May)
    PR.<u, x, y> = PolynomialRing(ZZ)
    Q = PR.quotient(x*y + 1 - u) # u = xy + 1
    polZ = Q(pol).lift()
 
    UU = XX*YY + 1
 
    # x-shifts
    gg = []
    for kk in range(mm + 1):
        for ii in range(mm - kk + 1):
            xshift = x**ii * modulus**(mm - kk) * polZ(u, x, y)**kk
            gg.append(xshift)
    gg.sort()
 
    # x-shifts list of monomials
    monomials = []
    for polynomial in gg:
        for monomial in polynomial.monomials():
            if monomial not in monomials:
                monomials.append(monomial)
    monomials.sort()
    
    # y-shifts (selected by Herrman and May)
    for jj in range(1, tt + 1):
        for kk in range(floor(mm/tt) * jj, mm + 1):
            yshift = y**jj * polZ(u, x, y)**kk * modulus**(mm - kk)
            yshift = Q(yshift).lift()
            gg.append(yshift) # substitution
    
    # y-shifts list of monomials
    for jj in range(1, tt + 1):
        for kk in range(floor(mm/tt) * jj, mm + 1):
            monomials.append(u**kk * y**jj)
 
    # construct lattice B
    nn = len(monomials)
    BB = Matrix(ZZ, nn)
    for ii in range(nn):
        BB[ii, 0] = gg[ii](0, 0, 0)
        for jj in range(1, ii + 1):
            if monomials[jj] in gg[ii].monomials():
                BB[ii, jj] = gg[ii].monomial_coefficient(monomials[jj]) * monomials[jj](UU,XX,YY)
 
    # Prototype to reduce the lattice
    if helpful_only:
        # automatically remove
        BB = remove_unhelpful(BB, monomials, modulus**mm, nn-1)
        # reset dimension
        nn = BB.dimensions()[0]
        if nn == 0:
            print("failure")
            return 0,0
 
    # check if vectors are helpful
    if debug:
        helpful_vectors(BB, modulus**mm)
    
    # check if determinant is correctly bounded
    det = BB.det()
    bound = modulus**(mm*nn)
    if det >= bound:
        print("We do not have det < bound. Solutions might not be found.")
        print("Try with highers m and t.")
        if debug:
            diff = (log(det) - log(bound)) / log(2)
            print("size det(L) - size e**(m*n) = ", floor(diff))
        if strict:
            return -1, -1
    else:
        print("det(L) < e**(m*n) (good! If a solution exists < N**delta, it will be found)")
 
    # display the lattice basis
    if debug:
        matrix_overview(BB, modulus**mm)
 
    # LLL
    if debug:
        print("optimizing basis of the lattice via LLL, this can take a long time")
 
    BB = reduce_lattice(BB, 0.8)
   # BB = BB.BKZ(block_size = 25)
    if debug:
        print("LLL is done!")
 
    # transform vector i & j -> polynomials 1 & 2
    if debug:
        print("looking for independent vectors in the lattice")
    found_polynomials = False
    
    for pol1_idx in range(nn - 1):
        for pol2_idx in range(pol1_idx + 1, nn):
            # for i and j, create the two polynomials
            PR.<w,z> = PolynomialRing(ZZ)
            pol1 = pol2 = 0
            for jj in range(nn):
                pol1 += monomials[jj](w*z+1,w,z) * BB[pol1_idx, jj] / monomials[jj](UU,XX,YY)
                pol2 += monomials[jj](w*z+1,w,z) * BB[pol2_idx, jj] / monomials[jj](UU,XX,YY)
 
            # resultant
            PR.<q> = PolynomialRing(ZZ)
            rr = pol1.resultant(pol2)

            # are these good polynomials?
            if rr.is_zero() or rr.monomials() == [1]:
                continue
            else:
                print("found them, using vectors", pol1_idx, "and", pol2_idx)
                found_polynomials = True
                break
        if found_polynomials:
            break
 
    if not found_polynomials:
        print("no independant vectors could be found. This should very rarely happen...")
        return 0, 0
    
    rr = rr(q, q)
    # solutions
    soly = rr.roots()
 
    if len(soly) == 0:
        print("Your prediction (delta) is too small")
        return 0, 0
 
    soly = soly[0][0]
    ss = pol1(q, soly)
    solx = ss.roots()[0][0]

    return solx, soly

def example1():
    # the modulus
    N = 0xa4d80845630d3b332f74f667ec8a0e49aba15b6f0c4f4006161d62c91b78cf6811421cc76609d2d9dba2c43be9d8ecdc6a0dff64a8041dcde52c7f92820b0a38fc91419e8ec9a5c69d47edc6e347934b4d87f97c5759886dac6c1143ff55b8eb11acfaa6cc70956a8ec7796e1a063b123bc2e467e30937c5a69c7ab5f8ed17e1
    # the public exponent
    e = 0x3458c2e97adef45f741c7db11ece6c0814aa5b6fad9144242cdaa16a6b4f3622477935f98a41765b92892b4de22a391cf08767447df113f5151c86edd109b97f9b045fd8ad5d7a51084684d4e2353db6c0e474d5d79f399a2bf4fd867ec85b7960845ab5497f705914912f797804c06dcff57139e040596d22b141e54835e0d3
 
    c = 0x91b097a5b1f6b12accdbda15cd2247384e1b3ed8311085a0f3e0dbb5fffce650a355600a02674189d1b7f4075df079c70354a08646e85ecf31dd150220cd1d4ce22d55a946500f4bd8def74fb0acea3e8d2e7bb1d27ebf2ca2e80fc28c3f0d88a041d4a556a18147f66b88c65f19c99b4b94c3f78d468b8accb4da7e7ce31b29
 
    # the hypothesis on the private exponent (the theoretical maximum is 0.292)
    delta = .250 # this means that d < N**delta
    print("赛题1")
    
    m = 4 # size of the lattice (bigger the better/slower)
    t = int((1-2*delta) * m)  # optimization from Herrmann and May
    X = 1 * floor(N**delta)  # this _might_ be too much
    Y = floor(N**(1/2))    # correct if p, q are ~ same size
    P.<x,y> = PolynomialRing(ZZ)
    A = int((N+1)/2)
    pol = 1 + x * (A + y)
 
    if debug:
        print("=== checking values ===")
        print("* delta:", delta)
        print("* delta < 0.292", delta < 0.292)
        print("* size of e:", int(log(e)/log(2)))
        print("* size of N:", int(log(N)/log(2)))
        print("* m:", m, ", t:", t)
 
    # boneh_durfee
    if debug:
        print("=== running algorithm ===")
        start_time = time.time()
 
    solx, soly = boneh_durfee(pol, e, m, t, X, Y)
    d = int(pol(solx, soly) / e)
    print("d is ", d)
    # found a solution?
    print("counting")
    if solx > 0:
        print("=== solution found ===")
        if False:
            print("x:", solx)
            print("y:", soly)
        if d > -2**268 and d < -2**267:
            print("identified here")
        #d = int(pol(solx, soly) / e)
        print("private key found:", d)
 
        # message = c.powermod(d, N)
        # print("message is ", message)
        # byte_string = bytes.fromhex(message.hex())
        # print("byte_string is ")
        # print(byte_string)
        # ascii_string = byte_string.decode("ASCII")
        # txt = ascii_string[::-1]
        # print(txt)
    else:
        print("=== no solution was found ===")
    if debug:
        print(("=== %s seconds ===" % (time.time() - start_time)))

def example2():
    # the modulus
    N = 0xd231f2c194d3971821984dec9cf1ef58d538975f189045ef8a706f6165aab4929096f61a3eb7dd8021bf3fdc41fe3b3b0e4ecc579b4b5e7e035ffcc383436c9656533949881dca67c26d0e770e4bf62a09718dbabc2b40f2938f16327e347f187485aa48b044432e82f5371c08f6e0bbde46c713859aec715e2a2ca66574f3eb
    # the public exponent
    e = 0x5b5961921a49e3089262761e89629ab6dff2da1504a0e5eba1bb7b20d63c785a013fd6d9e021c01baf1b23830954d488041b92bca2fe2c92e3373dedd7e625da11275f6f18ee4aef336d0637505545f70f805902ddbacb21bb8276d34a0f6dfe37ede87dd95bb1494dbb5763639ba3984240f1178e32aa36ee3c5fcc8115dde5
    c = 0x6a88a8fa2b8f28d96284298bab2061efeb35e3a086370e19523c15c429f5d783b9d4f32e31a402916f45ad4f2760ab30e77177335af44756bfbeef0f168b5e0dc8c3ddf75d141c358969cca0e7c2b8ab99ef8e33b031be1cbccd95b687682ac7b0dcc0d56f5651ee671d6358128d2e0801f247a6af4fe0dc5e8fb199eba0780f
    # the hypothesis on the private exponent (the theoretical maximum is 0.292)
    delta = .280 # this means that d < N**delta
    print("赛题 2")
    m = 13# size of the lattice (bigger the better/slower)
    t = int((1-2*delta) * m) # optimization from Herrmann and May
    X = 2 * floor(N**delta) # this _might_ be too much
    Y = floor(N**(1/2)) # correct if p, q are ~ same size
    P.<x,y> = PolynomialRing(ZZ)
    A = int((N+1)/2)
    pol = 1 + x * (A + y)
    if debug:
        print("=== checking values ===")
        print("* delta:", delta)
        print("* delta < 0.292", delta < 0.292)
        print("* size of e:", int(log(e)/log(2)))
        print("* size of N:", int(log(N)/log(2)))
        print("* m:", m, ", t:", t)
    # boneh_durfee
    if debug:
        print("=== running algorithm ===")
        start_time = time.time()
    solx, soly = boneh_durfee(pol, e, m, t, X, Y)
    d = int(pol(solx, soly) / e)
    print("d is ", d)
    # found a solution?
    print("counting")
    if solx > 0:
        print("=== solution found ===")
        if False:
            print("x:", solx)
            print("y:", soly)
        if d > -2**268 and d < -2**267:
            print("identified here")
        #d = int(pol(solx, soly) / e)
        print("private key found:", d)
        # message = c.powermod(d, N)
        # print("message is ", message)
        # byte_string = bytes.fromhex(message.hex())
        # print("byte_string is ")
        # print(byte_string)
        # ascii_string = byte_string.decode("ASCII")
        # txt = ascii_string[::-1]
        # print(txt)
    else:
        print("=== no solution was found ===")
    if debug:
        print(("=== %s seconds ===" % (time.time() - start_time)))


@parallel(ncpus=8) #537922560
def example3(msb):
    
    ############################################
    # How To Use This Script
    ##########################################

    #
    # The problem to solve (edit the following values)
    #r
    N = 0xf4c548636db62ffcc7ac4a0797952bea9a65bd426175af2435f72657e67ec8194667bfa94ce23c6f1e5baf3201867ab41701f6b8768e71009c41a3d5e9e7c109455341d549c7611f9f52851a2f017906aa9ccbedb95d238468e2c8577d30ecc4f158e3811fd5e2a6051443d468e3506bbc39bba710e34a604ac9e85d0feef8b3
    e = 0x16f4b438ba14e05afa944f7da9904f8c78ea52e4ca0be7fa2b5f84e22ddd7b0578a3477b19b7bb4a7f825acc45da2dd10e62dbd94a3386b97d92ee817b0c66c1507514a7860b9139bc2ac3a4e0fe304199214da00a4ca82bfcb7b18253e7e6144828e584dac2dfb9a03fabaf2376ce7c269923fbb60fc68325b9f6443e1f896f
    
    # the hypothesis on the private exponent (the theoretical maximum is 0.292)
    
    #print("赛题3 test")
    delta = 0.293
    m = 7
    s = 30 # number of msb enumerated
    t = ceil((1-2*delta) * m)
    pmsb = msb
    #print("pmsb is ", pmsb)
    qmsb = N/(pmsb * 2^(512 - s) + 2^(511 - s) )
    qmsb = floor(qmsb / 2^(512 - s))

    X = 2^300  # this _might_ be too much
    #Y = floor(N^(1/2))    # correct if p, q are ~ same size
    Y = floor(2*floor(N^(1/2)/2^s))
    #Y = 2^300
    #
    # Don't touch anything below
    #

    # Problem put in equation
    P.<x,y> = PolynomialRing(ZZ)
    #A = int((N+1)/2)
    A = N + 1-pmsb*2^(512 - s)-qmsb*2^(512 - s)
    pol = 1 + x * (A + y)

    #
    # Find the solutions!
    #

    # Checking bounds
    if debug:
        print("=== checking values ===")
        print("* delta:", delta)
        print("* delta < 0.292", delta < 0.292)
        print("* size of e:", int(log(e)/log(2)))
        print("* size of N:", int(log(N)/log(2)))
        print("* m:", m, ", t:", t)

    # boneh_durfee
    if debug:
        print("=== running algorithm ===")
        start_time = time.time()

    solx, soly = boneh_durfee(pol, e, m, t, X, Y)
    d0 = int(pol(solx, soly) / e)
    if solx > 0:
        return d0
    return d0

def test_LZQ23():
    # res = example(2^29)
    # print(res)
    # res = example([])

    start_time = time.time()
    flag = False
    #cnt = 0
    N = 0xf4c548636db62ffcc7ac4a0797952bea9a65bd426175af2435f72657e67ec8194667bfa94ce23c6f1e5baf3201867ab41701f6b8768e71009c41a3d5e9e7c109455341d549c7611f9f52851a2f017906aa9ccbedb95d238468e2c8577d30ecc4f158e3811fd5e2a6051443d468e3506bbc39bba710e34a604ac9e85d0feef8b3
    e = 0x16f4b438ba14e05afa944f7da9904f8c78ea52e4ca0be7fa2b5f84e22ddd7b0578a3477b19b7bb4a7f825acc45da2dd10e62dbd94a3386b97d92ee817b0c66c1507514a7860b9139bc2ac3a4e0fe304199214da00a4ca82bfcb7b18253e7e6144828e584dac2dfb9a03fabaf2376ce7c269923fbb60fc68325b9f6443e1f896f
    c = 0x26b1823cf836b226e2f5c90fdcd8420dbfcd02765b26e52ef3e5c0ab494c2f4650e475e280b0b5fff0d5016621186420b09e4706a5866e4a3319f23ef09d92c4e36acba39a0f6213fbe5ee1a736ce383e6e12351e6cbfd43f10a96b7fe34bdbaf948f2fb075d9063723c9f747fe6247ae9209e5d417faf2e37e6fee2eb863556
    g = 10 # Jumping gap
    s = 30 # number of msb enumerated
    parallel_num = 1
    for i in range(0, 2^g):
        print("i is ", i)
        if flag == True:
            break
        idx = 0
        while idx < 2^(s-g-1):
            print("idx is ", idx)
            if flag == True:
                break
            res = example3([(2^(s-1) + t * (2^g) + i) for t in range(idx, idx + parallel_num)]) #First bit is 1, so beginning at 2^(s-1)
            print("coming here")
            idx = idx + parallel_num
            for xx in res:
                #print(xx[1])
                if xx[1]!= 0:
                    print(xx)
                    flag = True
                    d = xx[1]
                    print("private key found:", d)
                    # message = c.powermod(d, N)
                    # print("message is ", message)
                    # byte_string = bytes.fromhex(message.hex())
                    # print("byte_string is ")
                    # print(byte_string)
                    # ascii_string = byte_string.decode("ASCII")
                    # txt = ascii_string[::-1]
                    break
                #print(xx[1])
            print(time.time() - start_time, " seconds")


if __name__ == "__main__":
    test_LZQ23()