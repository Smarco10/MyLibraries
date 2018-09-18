//
//  AESCrypto.h
//  test
//
//  Created by Marc-Antoine MARTIN on 15/06/2015.
//  Copyright (c) 2015 Marc-Antoine MARTIN. All rights reserved.
//

#ifndef AES_CRYPTO_H
#define AES_CRYPTO_H

#include <openssl/pem.h>
#include <openssl/ssl.h>
#include <openssl/evp.h>
#include <openssl/err.h>

//site: https://wiki.openssl.org/index.php/EVP_Symmetric_Encryption_and_Decryption

class AESCrypto{
private:
	unsigned char _key[SHA384_DIGEST_LENGTH];
    unsigned char _iv[EVP_MAX_IV_LENGTH];
	
#define AES_CYPHER_SIZE 16
	
	enum Modes{
		NONE,
		ENCRYPT,
		DECRYPT
	};
	
	Modes _mode;
	
	EVP_CIPHER_CTX *_ctx;
	unsigned long _outSize;
	int _curlen;
	unsigned char *_to;
	
public:
	AESCrypto(const unsigned char *key_p, size_t keylen_p);
	~AESCrypto(){
		/* Clean up */
		EVP_cleanup();
		memset(_key, 0, SHA384_DIGEST_LENGTH);
        memset(_iv, 0, EVP_MAX_IV_LENGTH);
	}
	
	void encrypt(unsigned char *to);
	void decrypt(unsigned char *to);

	unsigned long update(const unsigned char *from, unsigned int from_len);
	unsigned long finalize();
	void reset();
	
	inline unsigned long encrypt(const unsigned char *from, unsigned int fromlen, unsigned char *to){
		encrypt(to);
		
		if(!update(from, fromlen))
			return 0;
		
		return finalize();
	}
	
	inline unsigned long decrypt(const unsigned char *from, unsigned int fromlen, unsigned char *to){
		decrypt(to);
		
		if(!update(from, fromlen))
			return 0;
		
		return finalize();
	}
	
	static inline unsigned long encrypt(const unsigned char *key, unsigned int keylen, const unsigned char *from, unsigned int fromlen, unsigned char *to){
		if(key == NULL)
			return 0;
		return AESCrypto(key, keylen).encrypt(from, fromlen, to);
	}
	
	static inline unsigned long decrypt(const unsigned char *key, unsigned int keylen, const unsigned char *from, unsigned int fromlen, unsigned char *to){
		if(key == NULL)
			return 0;
		return AESCrypto(key, keylen).decrypt(from, fromlen, to);
	}
	
	static inline size_t getCypherSize(const size_t& inputSize){
		return inputSize + (AES_CYPHER_SIZE - (inputSize % AES_CYPHER_SIZE));
	}
	
	static inline void printError() {
        //ERR_print_errors_fp(stderr);
	}
};

#endif //AES_CRYPTO_H
