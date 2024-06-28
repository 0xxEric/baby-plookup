# A simple polynomial commitment interface

import random
import hashlib

from fft import fft
from ec import (G1Generator, G2Generator, default_ec)
from poly_utils import PrimeField
from pairing import ate_pairing

class PolyCommitment:
    def __init__(self, modulus,setup=None):
        assert setup == None
        self.modulus = modulus#默认值为default_ec.n，此处暂改为101，方便学习查看
        self.pf = PrimeField(self.modulus)  # order of Elliptic curve
        if setup == None:
            self.secret = random.randint(0, self.pf.modulus - 1)
            self.G1 = G1Generator()
            self.G2 = G2Generator()
            self.setup_vec1 = []
            self.setup_vec2 = []
        else:
            assert False


    def getSetupVector1(self, length):
        for i in range(len(self.setup_vec1), length):
            self.setup_vec1.append(self.G1 * (self.secret ** i)) #设χ=secert，则此处为1，χ，χ^2......χ^(n-1)
        return self.setup_vec1[0:length]

    def getSetupVector2(self, length):
        for i in range(len(self.setup_vec2), length):
            self.setup_vec2.append(self.G2 * (self.secret ** i)) #这里对于G2,也作了同样处理？与plonk文章中略有不同
        return self.setup_vec2[0:length]

    # get the commitment of a polynomial in evaluation form
    # return a curve point
    def getCommitment(self, evals, g):
        # g - roots of unity
        coeffs = fft(evals, self.pf.modulus, g, inv=True)
        sv = self.getSetupVector1(len(coeffs))
        return sum(s * c for s, c in zip(sv, coeffs))
    
    def getCommitmentbycoeffs(self, coeffs, g):
        sv = self.getSetupVector1(len(coeffs))
        return sum(s * c for s, c in zip(sv, coeffs))

    def getSingleProofByEvalIdx(self, evals, g, idx):
        # g - primitive root of unity
        # order - order of g
        # n - g^n root of unity
        # qx = (x^n - 1) / (x-x0)
        coeffs = fft(evals, self.pf.modulus, g, inv=True)
        y0 = evals[idx]
        x0 = self.pf.exp(g, idx)
        return self.getSingleProofAt(coeffs, x0, y0) #与getSingleProofAt功能类似，只是getSingleProofAt是在点（x0,y0)处打开，而本函数，是在第idx个根处打开，这里只是多一步把x0 y0对应得到

    def getSingleProofAt(self, coeffs, x0, y0):
        # g - primitive root of unity
        # order - order of g
        # n - g^n root of unity
        # qx = (x^n - 1) / (x-x0)
        coeffs = coeffs[:]
        coeffs[0] = coeffs[0] - y0
        qx = self.pf.div_polys(coeffs, [self.modulus-x0, 1]) #这里应该得到是商多项式qx的多项式(系数形式)，目测应该是：q(x)=(f(x)-y0)/(x-x0)，注意这里的q(x)应仍是系数形式
        sv = self.getSetupVector1(len(qx))
        return sum(s * c for s, c in zip(sv, qx)) #最后，得到的是多项式q(x)的承诺。可以认为即为多项式coeffs，在y0处打开为x0的承诺Cq(x)

    def verifySingleProof(self, commit, proof, x0, y0):
        sz2 = self.G2 * self.secret + (self.G2 * x0).negate() #即为：χ-ζ，ζ即为打开点x0
        pair0 = ate_pairing(proof, sz2) #这里proof，相当于是传入的Cq(x)

        cy = commit + (self.G1 * y0).negate() #这里cy即为：Cf(x)-y，y即对应打开点的y0
        pair1 = ate_pairing(cy, self.G2)
        return pair0 == pair1 #这里是验证单点打开


def test_poly_commitment(modulus):
    pc = PolyCommitment(modulus)
    G = pc.pf.exp(7, (pc.pf.modulus-1) // 4)  #待复核：这个G到底是什么？和g的关系是什么？ 落实：G=7^((modulus-1)/4)，用这个做生成元。类比于如果moduls=17，则使用G=（7^4)mod17=4做生成元，得到w0=1 w1=4 w2=16 w3=13,此时阶为（17-1）/4=4
    evals = [235, 2346, 132213, 61232]
    commit = pc.getCommitment(evals, G)
    proof = pc.getSingleProofByEvalIdx(evals, G, 1) # w^1，在x0=w处的打开
    assert pc.verifySingleProof(commit, proof, G, 2346) #验证，y0是不是等于f(x0)
    print("poly_commitment test passed")


def test_full_poly(modulus):
    # Verify full data blobs in EIP-4844 setup
    pc = PolyCommitment(modulus)
    G = pc.pf.exp(7, (pc.pf.modulus-1) // 4)
    evals = [235, 2346, 132213, 61232]
    commit = pc.getCommitment(evals, G)

    # given evals, check evals matches the commitment with single open
    # note that we could also do the check by
    # - find the evals coeffs
    # - multiple setup s^x G's
    # which may be expensive in contract.
    def get_proof(pc, commit, evals):
        # find a random evaluation point using Fiat-Shamir heuristic
        # where inputs are commit + evals
        data = bytes(commit)
        for eval in evals:
            # BLS modulus is in 256-bit
            data += eval.to_bytes(32, byteorder="big")  #把commit与evals拼接起来，即data
        # simple hash to point
        r = int.from_bytes(hashlib.sha256(data).digest(), byteorder="big") % default_ec.n #将拼接的data进行Hash，并取模，得到随机评估点r

        # use barycentric formula to calculate the point with evaluations
        pf = pc.pf
        yr = pf.eval_barycentric(r, [pf.exp(G, i) for i in range(4)], evals) #使用重心公式（Barycentric Formula）进行插值计算，计算在r处的取值

        # get the proof at random r (this part is off-chain)
        return pc.getSingleProofAt(fft(evals, pf.modulus, G, inv=True), r, yr) #然后，对随机点r的取值，进行验证（这里就与前面一样了）

    def verify_proof(pc, commit, evals, proof):
        # find a random evaluation point using Fiat-Shamir heuristic
        # where inputs are commit + evals
        data = bytes(commit)
        for eval in evals:
            # BLS modulus is in 256-bit
            data += eval.to_bytes(32, byteorder="big")
        # simple hash to point
        r = int.from_bytes(hashlib.sha256(data).digest(), byteorder="big") % default_ec.n

        # use barycentric formula to calculate the point with evaluations
        pf = pc.pf
        yr = pf.eval_barycentric(r, [pf.exp(G, i) for i in range(4)], evals)

        return pc.verifySingleProof(commit, proof, r, yr)

    proof = get_proof(pc, commit, evals)

    # single open to verify
    assert verify_proof(pc, commit, evals, proof)
    print("full_poly test passed")


if __name__ == "__main__":
    modulus=101
    test_poly_commitment(modulus)  #这里是只在一个点验证
    test_full_poly(modulus)
        
        
        
    