//
//  RSACrypto.h
//  test
//
//  Created by Marc-Antoine MARTIN on 17/05/2015.
//  Copyright (c) 2015 Marc-Antoine MARTIN. All rights reserved.
//

#ifndef RSA_CRYPTO_H
#define RSA_CRYPTO_H

#include <openssl/pem.h>
#include <openssl/rsa.h>
#include <openssl/ssl.h>
#include <openssl/evp.h>
#include <openssl/bio.h>
#include <openssl/err.h>

//site: http://hayageek.com/rsa-encryption-decryption-openssl-c/
//site cracker une clé RSA privée: http://www.parlonssecurite.com/factoriser-cracke-une-cle-rsa/

#define RSA_PUB_KEY_EXT "pub"
#define RSA_PRIV_KEY_EXT "priv"

typedef uint8_t RSAKeyInfo_t;
typedef int (err_callback_t)(const char*, size_t, void*);

class RSACrypto{
private:
	RSA *key = NULL;
    bool publicEncrypted = true;
    unsigned int rsa_max_size = 0;
    err_callback_t *err_cb = NULL;
    void *err_cb_u = NULL;
	
    #define RSA_PADDING_TYPE RSA_PKCS1_PADDING //RSA_PKCS1_OAEP_PADDING

    #if RSA_PADDING_TYPE == RSA_PKCS1_PADDING
		#define RSA_PADDING_SIZE RSA_PKCS1_PADDING_SIZE
    #else
		#define RSA_PADDING_SIZE 41
    #endif
	
    inline void printError() const {
        if(err_cb)
            ERR_print_errors_cb(err_cb, err_cb_u);
        else
            ERR_print_errors_fp(stderr);
	}
	
	static const RSAKeyInfo_t TYPE_BIT_IDLE = 0b00000000;
	
public:
	static const RSAKeyInfo_t TYPE_BIT_PRIV = 0b00000001;
	static const RSAKeyInfo_t TYPE_BIT_ENC = 0b00000010;
	static const RSAKeyInfo_t TYPE_BIT_VALID = 0b10000000;
	
	static RSAKeyInfo_t getKeyInfo(const char *keyFilename);
	static int getKeyInfo2(const char *keyFilename);

    static const unsigned int RSA_PRIV_KEY_DEFAULT_SIZE = 2048;
	
    RSACrypto(bool publicEncrypted = true);
	RSACrypto(const char *keyFilename, pem_password_cb *pass_cb = NULL, bool publicEncrypted = true);
    RSACrypto(unsigned int nbBits, bool publicEncrypted = true);
	~RSACrypto(){
		if(key != NULL)
			RSA_free(key);
		key = NULL;
	}
	
    bool setKey(const char* keyFileName, pem_password_cb *pass_cb);
    void setErrCallback(err_callback_t *err_cb, void *userdata){
        this->err_cb = err_cb;
        this->err_cb_u = userdata;
    }

    inline const RSA* getKey() const {return key;}
    inline bool canEncrypt() const { return key && ((publicEncrypted && key->n) || (!publicEncrypted && key->p && key->q)); }
    inline bool canDecrypt() const { return key && ((!publicEncrypted && key->n) || (publicEncrypted && key->p && key->q)); }
	
	inline unsigned int getEncryptMaxSize() const {
        return rsa_max_size - RSA_PADDING_SIZE;
    }

    inline unsigned int getDecryptMaxSize() const {
        return rsa_max_size;
    }
	
    int encrypt(const unsigned char *from, unsigned int fromLen, unsigned char *to) const;
    int decrypt(const unsigned char *from, unsigned int fromLen, unsigned char *to) const;
	
    void save(const char *publicKeyFilename, const char *privateKeyFilename, pem_password_cb *pass_cb = NULL) const;
};


#endif // RSA_CRYPTO_H
