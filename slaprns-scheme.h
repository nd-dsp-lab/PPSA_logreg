//
// Created by Antonia Januszewicz on 3/27/24.
//

#ifndef OPENFHE_SLAPRNS_H
#define OPENFHE_SLAPRNS_H

#include "PSA-base-scheme.h"

using namespace lbcrypto;

class SLAPScheme : public PSAScheme {
private:
    std::shared_ptr<CryptoParametersBFVRNS> cryptoParams;
    CCParams<CryptoContextCKKSRNS> CKKSparameters;
    CryptoContext<DCRTPoly> CKKSContext;

public:
    std::vector<BigInteger> delta_mod_q;
    std::vector<BigInteger> t_mod_q;
    std::vector<DoubleNativeInt> ComputeShit(DCRTPoly & plaintext);


    SLAPScheme(Scheme scheme, double scale);

    void Init();

    DCRTPoly Encrypt(const DCRTPoly& plaintext, const DCRTPoly &privateKey, const DCRTPoly &publicKey,
                             const bool do_noise,
                             double & noise_time, double & enc_time) override;

    DCRTPoly NSEncrypt(const DCRTPoly &plaintext, const DCRTPoly& privateKey, const DCRTPoly &publicKey);
    DCRTPoly MSEncrypt(const DCRTPoly& plaintext, const DCRTPoly& privateKey, const DCRTPoly& publicKey);

    DCRTPoly Decrypt(const std::vector<DCRTPoly>& ciphertexts, const DCRTPoly& aggregationKey, const uint64_t ts,
                     double & dec_time, unsigned int num_additions=0) override;

    DCRTPoly Decrypt(const std::vector<DCRTPoly>& ciphertexts, const DCRTPoly &aggregationKey, const DCRTPoly& publicKey,
                      double & dec_time, unsigned int num_additions=0) override;


    DCRTPoly NSDecrypt(const std::vector<DCRTPoly>& ciphertexts,const DCRTPoly& aggregationKey, const DCRTPoly &publicKey,
                       unsigned int num_additions=0);
    DCRTPoly MSDecrypt(const std::vector<DCRTPoly> &ciphertexts,const DCRTPoly& aggregationKey, const DCRTPoly &publicKey,
                       unsigned int num_additions=0);

    void SwitchBasis(DCRTPoly & ciphertext, DCRTPoly & plaintext);

    void ScaleDown(DCRTPoly & ciphertext, DCRTPoly & plaintext);

    DCRTPoly PolynomialEncrypt(const std::vector<double>& plaintext, const DCRTPoly &privateKey, const DCRTPoly& publicKey,
                               bool do_noise, double & noise_time,
                               double & enc_time, const std::vector<double>& e);

    std::vector<double> PolynomialDecrypt(const std::vector<DCRTPoly> &ciphertexts,std::vector<double> &constants, const DCRTPoly &aggregationKey, const DCRTPoly& publicKey,
                                          double & dec_time, unsigned int num_additions=0);

    std::vector<double> PolynomialDecrypt(const std::vector<DCRTPoly> &ciphertexts, std::vector<double> &constants, const DCRTPoly& aggregationKey, const uint64_t ts,
                                          double & dec_time, unsigned int num_additions=0);


    ~SLAPScheme() {};

};



#endif  //OPENFHE_SLAPRNS_H
