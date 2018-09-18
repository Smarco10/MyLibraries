//
//  RSACrypto.c
//
//  Created by Marc-Antoine MARTIN on 16/11/2015.
//  Copyright (c) 2015 Marc-Antoine MARTIN. All rights reserved.
//

#include "RSACrypto.h"
#include <fstream>

/**
 publicEncrypted = true -> encrypt with PUBLIC Key and decrypt with PRIVATE key (common usage)
 publicEncrypted = false -> encrypt with PRIVATE Key and decrypt with PUBLIC key (unsafe usage: discouraged)
 */
RSACrypto::RSACrypto(bool publicEncrypted):publicEncrypted(publicEncrypted){
	if(!(key = RSA_new()))
		return;
}

RSACrypto::RSACrypto(const char *keyFilename, pem_password_cb *pass_cb, bool publicEncrypted):RSACrypto(publicEncrypted){
    setKey(keyFilename, pass_cb);
}

RSACrypto::RSACrypto(unsigned int nbBits, bool publicEncrypted):publicEncrypted(publicEncrypted){
    BIGNUM *bne = NULL;
    unsigned long e = RSA_F4;

    // generate rsa key
    bne = BN_new();
    if(BN_set_word(bne, e) != 1){
        printError();
        BN_free(bne);
        bne = NULL;
        key = NULL;
        return;
    }

    key = RSA_new();
    if(RSA_generate_key_ex(key, nbBits, bne, NULL) != 1){
        printError();
        RSA_free(key);
        key = NULL;
        BN_free(bne);
        bne = NULL;
        return;
    }

    if(key->n != NULL)
        rsa_max_size = RSA_size(key);

    BN_free(bne);
    bne = NULL;
}

bool RSACrypto::setKey(const char* keyFilename, pem_password_cb *pass_cb){
    RSAKeyInfo_t kF = getKeyInfo(keyFilename);

    if(!(key && (kF & TYPE_BIT_VALID)))
        return false;

    BIO *bp = BIO_new_file(keyFilename, "r");

    if(!bp)
        return false;
	
	if(kF & TYPE_BIT_ENC && pass_cb)
		OpenSSL_add_all_ciphers();

    bool status = true;

    if(kF & TYPE_BIT_PRIV){
        //PRIVATE Key
        if(!PEM_read_bio_RSAPrivateKey(bp, &key, (kF & TYPE_BIT_ENC) ? pass_cb : NULL, (kF & TYPE_BIT_ENC) ? (void*)keyFilename : NULL)){
            printError();
            status = false;
        }
    } else {
        //PUBLIC Key
        if(!PEM_read_bio_RSA_PUBKEY(bp, &key, (kF & TYPE_BIT_ENC) ? pass_cb : NULL, (kF & TYPE_BIT_ENC) ? (void*)keyFilename : NULL)){
            printError();
            status = false;
        }
    }

    BIO_free_all(bp);

    if(key && key->n && status)
        rsa_max_size = RSA_size(key);

    return status;
}

RSAKeyInfo_t RSACrypto::getKeyInfo(const char *keyFilename){
	RSAKeyInfo_t type = TYPE_BIT_IDLE;
	
	char line[128] = {0};
	if(keyFilename != NULL){
		std::ifstream is(keyFilename);
		
		// KEY TYPE (PRIVATE | PUBLIC)
		is.getline(line, 128);
		
		if(strncmp(line, "-----BEGIN RSA PRIVATE KEY-----", 31) == 0)
			type |= TYPE_BIT_PRIV;
		else if(strncmp(line, "-----BEGIN PUBLIC KEY-----", 26) != 0){
			is.close();	
			return type;
		}
		
		// KEY (Proc-Type | DEK-Info | data)
		is.getline(line, 128);
		
		if(strncmp(line, "Proc-Type: 4,ENCRYPTED", 22) == 0)
			type |= TYPE_BIT_ENC;
		
		is.close();
		type |= TYPE_BIT_VALID;
	}
	
	return type;
}

int RSACrypto::encrypt(const unsigned char *from, unsigned int fromLen, unsigned char *to) const {
    if(key == NULL || !canEncrypt() || fromLen > getEncryptMaxSize() || from == NULL || to == NULL)
		return -1;
	
    int res = -1;

    if(publicEncrypted)
        res = RSA_public_encrypt(fromLen, from, to, key, RSA_PADDING_TYPE);
	else
        res = RSA_private_encrypt(fromLen, from, to, key, RSA_PADDING_TYPE);

    if(res == -1)
        printError();

    return res;
}

int RSACrypto::decrypt(const unsigned char *from, unsigned int fromLen, unsigned char *to) const {
    if(key == NULL || !canDecrypt() || fromLen > getDecryptMaxSize() || from == NULL || to == NULL)
		return -1;
	
    int res = -1;
	
    if(publicEncrypted)
        res = RSA_private_decrypt(fromLen, from, to, key, RSA_PADDING_TYPE);
	else
        res = RSA_public_decrypt(fromLen, from, to, key, RSA_PADDING_TYPE);

    if(res == -1)
        printError();

    return res;
}

void RSACrypto::save(const char *publicKeyFilename, const char *privateKeyFilename, pem_password_cb *pass_cb) const {
    if(!key)
        return;

    //save private key
    if(privateKeyFilename != NULL){
        BIO *bp_private = BIO_new_file(privateKeyFilename, "w+");
		PEM_write_bio_RSAPrivateKey(bp_private, key, pass_cb ? EVP_aes_256_cbc() : NULL, NULL, 0, pass_cb, pass_cb ? (void*)RSA_PRIV_KEY_EXT : NULL);
        BIO_free_all(bp_private);
    }

    //save public key
    if(publicKeyFilename != NULL){
        BIO *bp_public = BIO_new_file(publicKeyFilename, "w+");
        PEM_write_bio_RSA_PUBKEY(bp_public, key); //PEM_write_bio_RSAPublicKey
        BIO_free_all(bp_public);
	}
}
