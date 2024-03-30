//
// Created by Antonia Januszewicz on 3/27/24.
//

#ifndef OPENFHE_SLAPRNS_H
#define OPENFHE_SLAPRNS_H

#include "PSA-base-scheme.h"
#include <core/lattice/lat-hal.h>

using namespace lbcrypto;

class SLAPScheme : public PSAScheme {
public:
    DCRTPoly Encrypt(const DCRTPoly plaintext, const DCRTPoly privateKey,
                             bool do_noise,
                             double & noise_time, double & enc_time) override;

    DCRTPoly NSEncrypt(const DCRTPoly plaintext, const DCRTPoly privateKey, const DCRTPoly publicKey);
    DCRTPoly MSEncrypt(const DCRTPoly plaintext, const DCRTPoly privateKey, const DCRTPoly publicKey);

    DCRTPoly Decrypt(std::vector<DCRTPoly> ciphertexts, const DCRTPoly aggregationKey, const uint64_t ts,
                             double & dec_time, unsigned int num_additions=0) override;

    DCRTPoly NSDecrypt(const std::vector<DCRTPoly> ciphertexts,const DCRTPoly aggregationKey, const DCRTPoly publicKey,
                       unsigned int num_additions=0);
    DCRTPoly MSDecrypt(const std::vector<DCRTPoly> ciphertexts,const DCRTPoly aggregationKey, const DCRTPoly publicKey,
                       unsigned int num_additions=0);
};



#endif  //OPENFHE_SLAPRNS_H
