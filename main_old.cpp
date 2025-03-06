//#include <core/lattice/lat-hal.h>
//#include <pke/openfhe.h>
#include <iostream>
#include <getopt.h>
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
#include <math/dftransform.h>
#include "PSA-cryptocontext.h"

void handler(int sig) {
    void *array[10];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 10);

    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}

    int main(int argc, char ** argv) {
        signal(SIGSEGV, handler);
        std::cout << "Hello, World! " << std::endl;
        //DCRTPoly a = DCRTPoly();
        unsigned int plain_bits = 40; //log t
        unsigned int num_users = 2; //n
        unsigned int iters = 100; //i
        unsigned int k_prime = 1; //k
        Scheme scheme1 = NS;

        unsigned int N = 1; //N

        int c;
          while((c = getopt(argc, argv, "t:n:i:k:N:")) != -1){
            switch(c){
            case 't':{
                plain_bits = atoi(optarg);
                break;
            }
        case 'n':{
                num_users = atoi(optarg);
                break;
            }
        case 'i':{
                iters = atoi(optarg);
                break;
            }
        case 'k':{
                k_prime = atoi(optarg);
                break;
            }
        case 'N':{
                N = atoi(optarg);
                break;
            }
        default:{
            std::cout << "Invalid argument: " << c;
            if(optarg != nullptr){
                std::cout << ' ' << optarg;
            }
            std::cout << std::endl;
            return 1;
        }
            }
          }

        if(!plain_bits){
            throw std::runtime_error("Must have nonempty plaintext space");
        }  
        if(!num_users){
            throw std::runtime_error("Must have at least some users");
        }
        if(!iters){
            throw std::runtime_error("Must have at least some iterations");
        }

        unsigned int MAX_CTEXTS_DEFAULT = 20;

        //temp();

        //Code for testing SLAP, which isn't what this paper is about

        /**
        PSACryptocontext p = PSACryptocontext(plain_bits, num_users, iters, scheme1);
        std::vector<double> noise_times;
        std::vector<double> enc_times;
        std::vector<double> dec_times;
        p.TestEncryption(iters, false, noise_times, enc_times);

        p.TestDecryption(iters,dec_times);

        for(const double d : noise_times){
            std::cout << "noise_times " << d << '\n';
        }
        for(const double d : enc_times){
            std::cout << "enc_times " << d << '\n';
        }
        for(const double d : dec_times){
            std::cout << "dec_times " << d << '\n';
        }
         **/


        PSACryptocontext pp = PSACryptocontext(plain_bits, num_users, iters, scheme1);

        std::vector<double> poly_noise_times;
        std::vector<double> poly_enc_times;

        //pp.TestPolynomialEncryption(true, iters, poly_noise_times, poly_enc_times);
        // pp.TestPolynomialEncryption(1, MAX_CTEXTS_DEFAULT, poly_noise_times, poly_enc_times);

        pp.PolynomialEnvSetup(poly_noise_times, poly_enc_times);

        std::vector<double> inputvec(pp.aggregator.plaintextParams.GetRingDimension()/2,99);
        inputvec[2] = 1028;
        std::vector<double> expvec(pp.aggregator.plaintextParams.GetRingDimension()/2,1);

        std::cout << "First input: " << inputvec << std::endl;

        pp.PolynomialEncryption(inputvec, expvec, 1, poly_noise_times, poly_enc_times);

        std::vector<double> inputvec2(pp.aggregator.plaintextParams.GetRingDimension()/2,6);
        inputvec2[2] = 5;
        std::vector<double> expvec2(pp.aggregator.plaintextParams.GetRingDimension()/2,2);

        pp.PolynomialEncryption(inputvec, expvec2, 2, poly_noise_times, poly_enc_times);

        std::cout << "Second Input: " << inputvec2 << std::endl;

        // std::vector<double> inputvec3(pp.aggregator.plaintextParams.GetRingDimension()/2,6);
        // inputvec3[2] = 5;
        // std::vector<double> expvec3(pp.aggregator.plaintextParams.GetRingDimension()/2,2);
        // pp.PolynomialEncryption(inputvec3, expvec3, 3, poly_noise_times, poly_enc_times);

        // std::vector<double> inputvec4(pp.aggregator.plaintextParams.GetRingDimension()/2,6);
        // inputvec4[2] = 5;
        // std::vector<double> expvec4(pp.aggregator.plaintextParams.GetRingDimension()/2,2);
        // pp.PolynomialEncryption(inputvec4, expvec4, 4, poly_noise_times, poly_enc_times);

        // std::vector<double> inputvec5(pp.aggregator.plaintextParams.GetRingDimension()/2,6);
        // inputvec5[2] = 5;
        // std::vector<double> expvec5(pp.aggregator.plaintextParams.GetRingDimension()/2,2);
        // pp.PolynomialEncryption(inputvec5, expvec5, 5, poly_noise_times, poly_enc_times);

        // std::vector<double> inputvec6(pp.aggregator.plaintextParams.GetRingDimension()/2,6);
        // inputvec6[2] = 5;
        // std::vector<double> expve65(pp.aggregator.plaintextParams.GetRingDimension()/2,2);
        // pp.PolynomialEncryption(inputvec6, expvec6, 5, poly_noise_times, poly_enc_times);

        // std::vector<double> inputvec7(pp.aggregator.plaintextParams.GetRingDimension()/2,6);
        // inputvec7[2] = 5;
        // std::vector<double> expvec7(pp.aggregator.plaintextParams.GetRingDimension()/2,2);
        // pp.PolynomialEncryption(inputvec7, expvec7, 5, poly_noise_times, poly_enc_times);

        // std::vector<double> inputvec8(pp.aggregator.plaintextParams.GetRingDimension()/2,6);
        // inputvec8[2] = 5;
        // std::vector<double> expvec8(pp.aggregator.plaintextParams.GetRingDimension()/2,2);
        // pp.PolynomialEncryption(inputvec8, expvec8, 5, poly_noise_times, poly_enc_times);

        // std::vector<double> inputvec9(pp.aggregator.plaintextParams.GetRingDimension()/2,6);
        // inputvec9[2] = 5;
        // std::vector<double> expvec9(pp.aggregator.plaintextParams.GetRingDimension()/2,2);
        // pp.PolynomialEncryption(inputvec9, expvec9, 5, poly_noise_times, poly_enc_times);

        // std::vector<double> inputvec10(pp.aggregator.plaintextParams.GetRingDimension()/2,6);
        // inputvec10[2] = 5;
        // std::vector<double> expvec10(pp.aggregator.plaintextParams.GetRingDimension()/2,2);
        // pp.PolynomialEncryption(inputvec10, expvec10, 5, poly_noise_times, poly_enc_times);

        std::vector<double> decrypt_times;

        std::vector<double> constants(pp.aggregator.plaintextParams.GetRingDimension()/2,0.5);
        std::vector<double> outputvec = pp.PolynomialDecryption(constants, 1, decrypt_times);

        std::cout << "Final output: " << outputvec << std::endl;


        for(const double d : poly_noise_times){
            //std::cout << "poly_noise_times " << d << '\n';
        }
        for(const double d : poly_enc_times){
            //std::cout << "poly_enc_times " << d << '\n';
        }
        for(const double d : decrypt_times){
            //std::cout << "decrypt_times " << d << '\n';
        }

        return 0;
    }


