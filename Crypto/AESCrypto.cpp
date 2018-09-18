//
//  AESCrypto.c
//
//  Created by Marc-Antoine MARTIN on 16/11/2015.
//  Copyright (c) 2015 Marc-Antoine MARTIN. All rights reserved.

#include "AESCrypto.h"

AESCrypto::AESCrypto(const unsigned char *key_p, size_t keylen_p):_mode(NONE),_ctx(NULL),_outSize(0),_curlen(0),_to(NULL){
	
	/* Initialise the library */
    OpenSSL_add_all_algorithms();
	OpenSSL_add_all_ciphers();
	OpenSSL_add_all_digests();
	
	//Get 384-hash from key as: [256; 128] <=> [key, iv]
    //TODO: pb ne permet pas d'utiliser la ligne de commande pour décripter, voir EVP_BytesToKey()
    SHA384(key_p, keylen_p, _key);
    memcpy(_iv, _key + (SHA384_DIGEST_LENGTH - EVP_MAX_IV_LENGTH), EVP_MAX_IV_LENGTH);
}

void AESCrypto::encrypt(unsigned char *to) {
	//clean all
	reset();

	/* Create and initialise the context */
	if(!(_ctx = EVP_CIPHER_CTX_new())){
		printError();
		return;
	}
	
	/* Initialise the encryption operation. IMPORTANT - ensure you use a key
	 * and IV size appropriate for your cipher
	 * In this example we are using 256 bit AES (i.e. a 256 bit key). The
	 * IV size for *most* modes is the same as the block size. For AES this
	 * is 128 bits */
    if(1 != EVP_EncryptInit_ex(_ctx, EVP_aes_256_cbc(), NULL, _key, _iv)){
		printError();
		return;
	}
	
	_to = to;
	_mode = ENCRYPT;
}

void AESCrypto::decrypt(unsigned char *to) {
	
	//clean all
	reset();
	
  /* Create and initialise the context */
	if(!(_ctx = EVP_CIPHER_CTX_new())){
		printError();
		return;
	}
	
  /* Initialise the decryption operation. IMPORTANT - ensure you use a key
   * and IV size appropriate for your cipher
   * In this example we are using 256 bit AES (i.e. a 256 bit key). The
   * IV size for *most* modes is the same as the block size. For AES this
   * is 128 bits */
	if(1 != EVP_DecryptInit_ex(_ctx, EVP_aes_256_cbc(), NULL, _key, _iv)){
		printError();
		return;
	}
	
	_to = to;
	_mode = DECRYPT;
}

unsigned long AESCrypto::update(const unsigned char *from, unsigned int from_len){
	
	if(!_ctx || !_to || !from)
		return 0;

    unsigned char *fromB;
    unsigned int from_lenB = from_len;

    if(from_len < AES_CYPHER_SIZE){
        fromB = new unsigned char[AES_CYPHER_SIZE];
        memcpy(fromB, from, from_len);
        memset(fromB + from_len, 0, AES_CYPHER_SIZE - from_len);
        from_lenB = AES_CYPHER_SIZE;
    }

    bool failed = false;
	
	switch(_mode){
		case ENCRYPT:
            if(1 != EVP_EncryptUpdate(_ctx, _to, &_curlen, from_len < AES_CYPHER_SIZE ? fromB : from, from_lenB)){
				printError();
                failed = true;
			}
			
			break;
		case DECRYPT:
            if(1 != EVP_DecryptUpdate(_ctx, _to, &_curlen, from_len < AES_CYPHER_SIZE ? fromB : from, from_lenB)){
				printError();
                failed = true;
			}
			break;
		case NONE:
		default:
            failed = true;
			break;
	}

    if(from_len < AES_CYPHER_SIZE){
        delete[] fromB;
    }
	
    if(failed)
        return 0;

	_outSize += _curlen;
	return _outSize;
	
}

unsigned long AESCrypto::finalize(){
	
	/* Finalise the decryption. Further plaintext bytes may be written at
	 * this stage.
	 */
	
	if(!_ctx || !_to)
		return 0;
	
	switch(_mode){
		case ENCRYPT:
			if(1 != EVP_EncryptFinal_ex(_ctx, _to + _curlen, &_curlen)){
				printError();
				return 0;
			}
			break;
		case DECRYPT:
			if(1 != EVP_DecryptFinal_ex(_ctx, _to + _curlen, &_curlen)){
				printError();
				return 0;
			}
			break;
		case NONE:
		default:
			return 0;
			break;
	}
	
	_outSize += _curlen;
	unsigned long outSize = _outSize;
	reset();
	
	return outSize;
}

void AESCrypto::reset(){
	_mode = NONE;
	
	/* Clean up */
	if(_ctx != NULL)
		EVP_CIPHER_CTX_free(_ctx);
	
	_ctx = NULL;
	_outSize = 0;
	_curlen = 0;
	_to = NULL;
}
