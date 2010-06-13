/*
   +----------------------------------------------------------------------+
   | This program is free software; you can redistribute it and/or        |
   | modify it under the terms of the GNU General Public License          |
   | as published by the Free Software Foundation; either version 3       |                 
   | of the License, or (at your option) any later version.               | 
   |                                                                      |
   | This program is distributed in the hope that it will be useful,      |
   | but WITHOUT ANY WARRANTY; without even the implied warranty of       |
   | MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        |
   | GNU General Public License for more details.                         | 
   |                                                                      |
   | You should have received a copy of the GNU General General Public    |
   | License in the file LICENSE along with this library;                 |
   | if not, write to the                                                 | 
   |                                                                      |
   |   Free Software Foundation, Inc.,                                    |
   |   51 Franklin Street, Fifth Floor,                                   |
   |   Boston, MA  02110-1301  USA                                        |
   +----------------------------------------------------------------------+
   | Authors: Evan Klinger <eklinger@phash.org>                           |
   +----------------------------------------------------------------------+
*/

/* $ Id: $ */ 

#ifndef PHP_PHASH_H
#define PHP_PHASH_H

#ifdef  __cplusplus
extern "C" {
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <php.h>

#ifdef HAVE_PHASH

#include <php_ini.h>
#include <SAPI.h>
#include <ext/standard/info.h>
#include <Zend/zend_extensions.h>
#ifdef  __cplusplus
} // extern "C" 
#endif
#include <pHash.h>
#include <audiophash.h>
#ifdef  __cplusplus
extern "C" {
#endif

extern zend_module_entry pHash_module_entry;
#define phpext_pHash_ptr &pHash_module_entry

#ifdef PHP_WIN32
#define PHP_PHASH_API __declspec(dllexport)
#else
#define PHP_PHASH_API
#endif

PHP_MINIT_FUNCTION(pHash);
PHP_MSHUTDOWN_FUNCTION(pHash);
PHP_RINIT_FUNCTION(pHash);
PHP_RSHUTDOWN_FUNCTION(pHash);
PHP_MINFO_FUNCTION(pHash);

#ifdef ZTS
#include "TSRM.h"
#endif

#define FREE_RESOURCE(resource) zend_list_delete(Z_LVAL_P(resource))

#define PROP_GET_LONG(name)    Z_LVAL_P(zend_read_property(_this_ce, _this_zval, #name, strlen(#name), 1 TSRMLS_CC))
#define PROP_SET_LONG(name, l) zend_update_property_long(_this_ce, _this_zval, #name, strlen(#name), l TSRMLS_CC)

#define PROP_GET_DOUBLE(name)    Z_DVAL_P(zend_read_property(_this_ce, _this_zval, #name, strlen(#name), 1 TSRMLS_CC))
#define PROP_SET_DOUBLE(name, d) zend_update_property_double(_this_ce, _this_zval, #name, strlen(#name), d TSRMLS_CC)

#define PROP_GET_STRING(name)    Z_STRVAL_P(zend_read_property(_this_ce, _this_zval, #name, strlen(#name), 1 TSRMLS_CC))
#define PROP_GET_STRLEN(name)    Z_STRLEN_P(zend_read_property(_this_ce, _this_zval, #name, strlen(#name), 1 TSRMLS_CC))
#define PROP_SET_STRING(name, s) zend_update_property_string(_this_ce, _this_zval, #name, strlen(#name), s TSRMLS_CC)
#define PROP_SET_STRINGL(name, s, l) zend_update_property_stringl(_this_ce, _this_zval, #name, strlen(#name), s, l TSRMLS_CC)


#if HAVE_VIDEO_HASH
PHP_FUNCTION(ph_dct_videohash);
#if (PHP_MAJOR_VERSION >= 5)
ZEND_BEGIN_ARG_INFO_EX(ph_dct_videohash_arg_info, ZEND_SEND_BY_VAL, ZEND_RETURN_VALUE, 1)
  ZEND_ARG_INFO(0, file)
ZEND_END_ARG_INFO()
#else /* PHP 4.x */
#define ph_dct_videohash_arg_info NULL
#endif

#endif /* HAVE_VIDEO_HASH */
#if HAVE_IMAGE_HASH
PHP_FUNCTION(ph_dct_imagehash);
#if (PHP_MAJOR_VERSION >= 5)
ZEND_BEGIN_ARG_INFO_EX(ph_dct_imagehash_arg_info, ZEND_SEND_BY_VAL, ZEND_RETURN_VALUE, 1)
  ZEND_ARG_INFO(0, file)
ZEND_END_ARG_INFO()
#else /* PHP 4.x */
#define ph_dct_imagehash_arg_info NULL
#endif

#endif /* HAVE_IMAGE_HASH */
PHP_FUNCTION(ph_texthash);
#if (PHP_MAJOR_VERSION >= 5)
ZEND_BEGIN_ARG_INFO_EX(ph_texthash_arg_info, ZEND_SEND_BY_VAL, ZEND_RETURN_VALUE, 1)
  ZEND_ARG_INFO(0, file)
ZEND_END_ARG_INFO()
#else /* PHP 4.x */
#define ph_texthash_arg_info NULL
#endif

#if HAVE_AUDIO_HASH
PHP_FUNCTION(ph_audiohash);
#if (PHP_MAJOR_VERSION >= 5)
ZEND_BEGIN_ARG_INFO_EX(ph_audiohash_arg_info, ZEND_SEND_BY_VAL, ZEND_RETURN_VALUE, 1)
  ZEND_ARG_INFO(0, file)
  ZEND_ARG_INFO(0, sample_rate)
  ZEND_ARG_INFO(0, channels)
ZEND_END_ARG_INFO()
#else /* PHP 4.x */
#define ph_audiohash_arg_info NULL
#endif

#endif /* HAVE_AUDIO_HASH */
#if HAVE_IMAGE_HASH
PHP_FUNCTION(ph_image_dist);
#if (PHP_MAJOR_VERSION >= 5)
ZEND_BEGIN_ARG_INFO_EX(ph_image_dist_arg_info, ZEND_SEND_BY_VAL, ZEND_RETURN_VALUE, 2)
  ZEND_ARG_INFO(0, h1)
  ZEND_ARG_INFO(0, h2)
ZEND_END_ARG_INFO()
#else /* PHP 4.x */
#define ph_image_dist_arg_info NULL
#endif

#endif /* HAVE_IMAGE_HASH */
#if HAVE_VIDEO_HASH
PHP_FUNCTION(ph_video_dist);
#if (PHP_MAJOR_VERSION >= 5)
ZEND_BEGIN_ARG_INFO_EX(ph_video_dist_arg_info, ZEND_SEND_BY_VAL, ZEND_RETURN_VALUE, 2)
  ZEND_ARG_INFO(0, h1)
  ZEND_ARG_INFO(0, h2)
  ZEND_ARG_INFO(0, thresh)
ZEND_END_ARG_INFO()
#else /* PHP 4.x */
#define ph_video_dist_arg_info NULL
#endif

#endif /* HAVE_VIDEO_HASH */
#if HAVE_AUDIO_HASH
PHP_FUNCTION(ph_audio_dist);
#if (PHP_MAJOR_VERSION >= 5)
ZEND_BEGIN_ARG_INFO_EX(ph_audio_dist_arg_info, ZEND_SEND_BY_VAL, ZEND_RETURN_VALUE, 2)
  ZEND_ARG_INFO(0, h1)
  ZEND_ARG_INFO(0, h2)
  ZEND_ARG_INFO(0, block_size)
  ZEND_ARG_INFO(0, thresh)
ZEND_END_ARG_INFO()
#else /* PHP 4.x */
#define ph_audio_dist_arg_info NULL
#endif

#endif /* HAVE_AUDIO_HASH */
PHP_FUNCTION(ph_compare_text_hashes);
#if (PHP_MAJOR_VERSION >= 5)
ZEND_BEGIN_ARG_INFO_EX(ph_compare_text_hashes_arg_info, ZEND_SEND_BY_VAL, ZEND_RETURN_VALUE, 2)
  ZEND_ARG_INFO(0, h1)
  ZEND_ARG_INFO(0, h2)
ZEND_END_ARG_INFO()
#else /* PHP 4.x */
#define ph_compare_text_hashes_arg_info NULL
#endif

#ifdef  __cplusplus
} // extern "C" 
#endif

#endif /* PHP_HAVE_PHASH */

#endif /* PHP_PHASH_H */


/*
 * Local variables:
 * tab-width: 4
 * c-basic-offset: 4
 * End:
 * vim600: noet sw=4 ts=4 fdm=marker
 * vim<600: noet sw=4 ts=4
 */
