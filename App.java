import java.math.BigInteger;
import java.util.Random;

public class App {
    //REQUIRED FUNCTION IN THE ASSIGNMENT

    /* Function 1: Find large prime number when given the number of bits
    Input: num_bits: number of bit
    Output: large prime number
    Algorithm:  1. Generate a random large number
                2. Check if that is a prime number (Using function: isPrime(number))
                3. If True: return large prime number found
                4. If False: Loop until finding a large prime number
    */ 
    static BigInteger getBigIntegerPrimeNumber(int num_bits) {
        Random rand = new Random();
        BigInteger primeCandidate;
        do {
            primeCandidate = new BigInteger(num_bits, rand);
        } while (!isPrime(primeCandidate));
        return primeCandidate;
    }
    /* Function 2: Calculate the greatest common divisor when given two large number
    Input: two large numbers: a and b
    Output: the greatest common divisor of a and b
    Algorithm:  Using Euclid Algorithm
    */
    static BigInteger getGCD(BigInteger a, BigInteger b) {
        while (true) {
            if (a.equals(BigInteger.ZERO)) return b;
            if (b.equals(BigInteger.ZERO)) return a;
            if (a.compareTo(b) == 1) {
                a = a.mod(b);
            } else {
                b = b.mod(a);
            }
        }
    }
    /* Function 3: Calculate the decryption key d given the encryption key e and two large prime numbers
    Input: encryption key e, large prime number p, q. (e and (p-1)*(q-1) are two coprime prime numbers)
    Output: decryption key d
    Algorithm:  1. Calculate the φ(n) = (p-1)*(q-1)
                2. Calculate d = Modular multiplicative inverse (φ(n)) --- Using Extended Euclidean Algorithm
    */
    static BigInteger findDecryptionKey_d(BigInteger e, BigInteger p, BigInteger q) {
        BigInteger phi_n = (p.subtract(BigInteger.ONE)).multiply(q.subtract(BigInteger.ONE));
        BigInteger d = mod_Inverse(e,phi_n);
        return d;
    }
    /* Function 4:     // Generate a random key pair given two large prime numbers
    Input: large prime number p, q. 
    Output: key pair (e,d) 
    Algorithm:  1. Calculate the φ(n) = (p-1)*(q-1)
                2. Loop and choose a large number e randomly such that e and φ(n) are two coprime integers.
                3. Calculate d = Modular multiplicative inverse of a number e modulo 𝜙(𝑛) --- Using Extended Euclidean Algorithm 
    */
    static BigInteger[] generateRandomKeyPair(BigInteger p, BigInteger q) {
        BigInteger phi_n = p.subtract(BigInteger.ONE).multiply(q.subtract(BigInteger.ONE));
        Random random = new Random();
        BigInteger e = new BigInteger(phi_n.bitLength(), random);
        while (e.compareTo(BigInteger.ONE) <= 0 || e.compareTo(phi_n) >= 0 || !e.gcd(phi_n).equals(BigInteger.ONE)) {
            e = new BigInteger(phi_n.bitLength(), random);
        }
        BigInteger d = mod_Inverse(e,phi_n);
        return new BigInteger[]{e, d};
    }
    /* Function 5:     // Encrypt the plain-text input and return the cipher-text
    Input: Public key of RSA: public large number e; large prime number n; plain-text m
    Output: Cipher-text c
    Algorithm:  Calculate c = m^e mod(n)
    */
    static BigInteger Encrypt(BigInteger m, BigInteger e, BigInteger n) {
        return modulo_Pow(m, e, n);
    }
    /* Function 6:     // Decrypt the cipher-text input and return the plain-text
    Input: Private key of RSA: private key d; large prime number n; cipher-text c
    Output: Plain-text m
    Algorithm:  Calculate m = c^d mod(n)
    */
    static BigInteger Decrypt(BigInteger c, BigInteger d, BigInteger n) {
        return modulo_Pow(c, d, n);
    }

//---------------------------------------------------------------------------------------------------------------------------------------
//HELPER FUNCTION 
        /*
        Function mod_Inverse:  Receive e, n -> return d such that e*d mod(n) = 1
        Algorithm:  1. Check to ensure that e and n are coprime integers.
                    2. Find d through Euclidean Extended Algorithm
                    3. Ensure that d>0 by add d with n until 0<d<n
        */

    static BigInteger mod_Inverse(BigInteger e, BigInteger n) {
        
    if (!getGCD(e, n).equals(BigInteger.ONE)) {
            throw new ArithmeticException("Please check to ensure that e and n are coprime integers");
        }
        BigInteger d = extendedEuclidean(e, n);
        while (d.compareTo(BigInteger.ZERO) < 0) {
            d = d.add(n);
        }
        return d;
    }

    /*
        Function: extendedEuclian receive large number a and b; Calculate: GCD(a,b),x,y
        Algorithm:  1. Algorithm will find x, y such that ax+by=GCD(a,b)
                    2. With input a,b called from mod_Inverse, a <- e, b <- n => a and b are co-prime => GCD(a,b)=1
                    3. We have ax+by = 1 => (ax+by) mod (b) = 1 => ax mod (b) = 1 => x is the Modular multiplicative inverse of number (a) modulo (b)
                    4. Return x
    */

    private static BigInteger extendedEuclidean(BigInteger a, BigInteger b) {
        BigInteger x = BigInteger.ZERO;
        BigInteger y = BigInteger.ONE;
        BigInteger lastX = BigInteger.ONE;
        BigInteger lastY = BigInteger.ZERO;
        BigInteger temp;
        while (!b.equals(BigInteger.ZERO)) {
            BigInteger quotient = a.divide(b);
            BigInteger remainder = a.mod(b);
            a = b;
            b = remainder;
            temp = x;
            x = lastX.subtract(quotient.multiply(x));
            lastX = temp;
            temp = y;
            y = lastY.subtract(quotient.multiply(y));
            lastY = temp;
        }
        return lastX;
    }
    /*
        Function: Return a^b mod(n) with a and b are 2 large number
        Algorithm:  1. result = 1
                    2. a = a mod (n)
                    3. Loop while b >0: 
                        3.1. If b is odd -> result = (result*a) mod (n)
                        3.2. b <- b/2 and a = (a*a) mod (n)
                    4. Return result
    */
    static BigInteger modulo_Pow(BigInteger a, BigInteger b, BigInteger n) {
        BigInteger result = BigInteger.ONE;
        a = a.mod(n);
        while (b.compareTo(BigInteger.ZERO) > 0) { 
            if (b.and(BigInteger.ONE).equals(BigInteger.ONE)) { 
                result = result.multiply(a).mod(n); 
            }
            b = b.shiftRight(1);
            a = a.multiply(a).mod(n); 
        }
        return result;
    }
    /*
        Function: IsPrime: Receive a large number n, return True if n is prime else return False
        Algorithm:  1. if n =1 or 2: n is prime; if n is even and n!=2: n is not a prime
                    2. Calculate n-1 to d*2^r
                    3. Checking Miller_Rabin Loop with k = 20 times 
                    4. If n pass 20 loop return True, else: False
    */
    static boolean isPrime(BigInteger n) {
        if (n.compareTo(BigInteger.ONE) <= 0) {
            return false; 
        }
        if (n.equals(BigInteger.TWO)) {
            return true; 
        }
        if (n.mod(BigInteger.TWO).equals(BigInteger.ZERO)) {
            return false;
        }

        BigInteger d = n.subtract(BigInteger.ONE);
        int r = 0;
        while (d.mod(BigInteger.TWO).equals(BigInteger.ZERO)) {
            d = d.divide(BigInteger.TWO);
            r++;
        }

        int k = 20;
        for (int i = 0; i < k; i++) {
            if (!Miller_Rabin(n, d, r)) {
                return false;
            }
        }
        return true;
    }

    /*
        Function: Miller_Rabin: Receive 2 large number n, d and int r; return True if n pass the test of Miller_Rabin Test
        Algorithm:  1. Choose number a randomly in range [2,n-2]
                    2. Calculate n-1 to d*2^r
                    3. x = a^d mod n, if x = 1 or x = n - 1 return True (follow the Miller_Rabin Algorithm - n MAY BE a prime)
                    4. For j=0 to r-1:
                            x <- x^2 mod n; if x = 1 or x = n - 1 return True (at here it is the same with a^(d*2*j) for j from 0 to r-1)
    */
    private static boolean Miller_Rabin(BigInteger n, BigInteger d, int r) {
        Random rand = new Random();
        BigInteger a = new BigInteger(n.bitLength(), rand);
        a = a.mod(n.subtract(BigInteger.TWO)).add(BigInteger.TWO); 
        BigInteger x = modulo_Pow(a,d, n);
        if (x.equals(BigInteger.ONE) || x.equals(n.subtract(BigInteger.ONE))) {
            return true;
        }
        for (int i = 0; i < r - 1; i++) {
            x = modulo_Pow(x,BigInteger.TWO, n);
            if (x.equals(BigInteger.ONE) || x.equals(n.subtract(BigInteger.ONE))) {
                return true;
            }
        }
        return false;
    }

    /*
        In cryptography, a prime number p is said to be "strong" if the following conditions are satisfied.
            1. p is sufficiently large to be useful in cryptography; typically this requires p to be too large 
                for plausible computational resources to enable a cryptanalyst to factorise products of p with other strong primes.
            2. p − 1 has large prime factors. That is, p = a1q1 + 1 for some integer a1 and large prime q1. (Prevent factoring attack Pollard's p − 1 algorithm)
            3. q1 − 1 has large prime factors. That is, q1 = a2q2 + 1 for some integer a2 and large prime q2. (Prevent reducing the odds that an RSA cycling attack succeeds)
            4. p + 1 has large prime factors. That is, p = a3q3 − 1 for some integer a3 and large prime q3. (Prevent factoring attack Williams's p + 1 algorithm)
    */



    /*
    Function: findStrongPrime (int bit); receive number of bit, return strong prime number.
    Algorithm: 
            1. generate a strong level 2 prime number q1 such that q1 = a2q2+1 (follow rule 3)
            2. generate a random prime number q3 
            3. Calculate and return strong prime number p such that p = a1q1 + 1 and p = a3q3 − 1 (a1 and a3 are integers)

    */
    static BigInteger findStrongPrime(int bit) {
        BigInteger two = new BigInteger("2");
        Random rand = new Random();
        BigInteger res = new BigInteger(bit,rand);
        BigInteger p0 = findStrongPrime_p1(bit/2);
        BigInteger p1 = getBigIntegerPrimeNumber(bit/2);
        while (p1.compareTo(p0) == 0) {
            p1 = getBigIntegerPrimeNumber(bit/2);
        }
        p1 = p1.multiply(two);

        BigInteger increment = p0.multiply(p1);

        BigInteger inv1 = mod_Inverse(p1, p0);
        BigInteger crt1 = inv1.multiply(p1);
        BigInteger inv2 = mod_Inverse(p0, p1);
        BigInteger crt2 = p1.subtract(inv2); 
        crt2 = crt2.multiply(p0);

        BigInteger crt = crt1.add(crt2).add(increment);
        BigInteger resmod = res.mod(increment);
        res = res.add(crt.subtract(resmod).mod(increment));
        increment = increment.multiply(two);

        while (!isPrime(res)) {
            res = res.add(increment);
        }
        return res;
     }
     /*
    Function: findStrongPrime (int bit); receive number of bit, return strong prime number level 2.
    Algorithm: 
            1. generate a random prime number q
            2. Calculate and return strong prime number level 2 p such that p = aq + 1 (follow rule 3)

    */
     static BigInteger findStrongPrime_p1(int bit) {
        BigInteger two = new BigInteger("2");

        Random rand = new Random();
        BigInteger res = new BigInteger(bit,rand);
        BigInteger p = getBigIntegerPrimeNumber(bit/2);
        p = p.multiply(two);
        BigInteger increment = p.multiply(p);

        BigInteger crt = mod_Inverse(BigInteger.ONE, p);
        crt = crt.add(increment);

        BigInteger resmod = res.mod(increment);
        res = res.add(crt.subtract(resmod).mod(increment));
        increment = increment.multiply(two);

        while (!isPrime(res)) {
            res = res.add(increment);
        }
        return res;
     }
    
//----------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------- RUN AND EVALUATE ------------------------------------------------------------------------

//  Run the RSA_Encryption using normal random prime number generator function.
    static void normalRSA_Encryption(int num_bits){ 
        // Create key
        BigInteger q;
        BigInteger p = getBigIntegerPrimeNumber(num_bits);
        do{
         q = getBigIntegerPrimeNumber(num_bits);}
        while (q.compareTo(p)==0);
        BigInteger [] pair = generateRandomKeyPair(p, q);
        BigInteger e = pair[0];
        BigInteger d = pair[1];
        BigInteger n = p.multiply(q);
        System.out.println("NORMAL RSA-ENCRYPTION");System.out.println("Public key");
        System.out.println("\te: "+e);
        System.out.println("\tn: "+n);
        System.out.println("Private key");
        System.out.println("\td: "+d);
        System.out.println("\tp: "+p);
        System.out.println("\tq: "+q);
        Random rand = new Random();
        BigInteger plain_text = new BigInteger(num_bits*2, rand);
        plain_text = plain_text.mod(p.multiply(q));
        System.out.println("Plain_text:");
        System.out.println(plain_text);
        System.out.println("Encryption:");
        BigInteger c = Encrypt(plain_text, e, n);
        System.out.println(c);  
        System.out.println("Decryption:");
        plain_text = Decrypt(c, d, n);
        assert(plain_text.compareTo(p.multiply(q)) == -1);
        System.out.println(plain_text);
    }
    //---------------------- IMPROVED RSA ENCRYPTION USING STRONG PRIME NUMBER----------------------------------------------------
    static void improvedRSA_Encryption(int num_bits){
        // Create key
        BigInteger q;
        BigInteger e;
        BigInteger d;
        BigInteger p = findStrongPrime(num_bits);
        do{
         q = findStrongPrime(num_bits);}
        while (q.compareTo(p)==0);
        do{
        BigInteger [] pair = generateRandomKeyPair(p, q);
         e = pair[0];
         d = pair[1];
        } while (!getGCD(e, d).equals(BigInteger.ONE));
        BigInteger n = p.multiply(q);
        System.out.println("IMPROVED RSA-ENCRYPTION");System.out.println("Public key");
        System.out.println("\te: "+e);
        System.out.println("\tn: "+n);
        System.out.println("Private key");
        System.out.println("\td: "+d);
        System.out.println("\tp: "+p);
        System.out.println("\tq: "+q);
        Random rand = new Random();
        BigInteger plain_text = new BigInteger(num_bits*2, rand);
        plain_text = plain_text.mod(p.multiply(q));
        System.out.println("Plain_text:");
        System.out.println(plain_text);
        System.out.println("Encryption:");
        BigInteger c = Encrypt(plain_text, e, n);
        System.out.println(c);  
        System.out.println("Decryption:");
        plain_text = Decrypt(c, d, n);
        assert(plain_text.compareTo(p.multiply(q)) == -1);
        System.out.println(plain_text);
    }
    

    public static void main(String args[]) throws Exception {
        long startTime = System.currentTimeMillis();
        // normalRSA_Encryption(2048);
        improvedRSA_Encryption(2048);
        long endTime = System.currentTimeMillis();
        long elapsedTime = endTime - startTime;
        System.out.println("Running time of the RSA-Encryption: " + elapsedTime + " milliseconds");
       
    }

}