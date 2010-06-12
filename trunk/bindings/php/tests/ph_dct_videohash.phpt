--TEST--
ph_dct_videohash() function
--SKIPIF--
<?php 

if(!extension_loaded('pHash')) die('skip ');

if(!function_exists('ph_dct_videohash')) die('skip not compiled in (HAVE_VIDEO_HASH)');

 ?>
--FILE--
<?php
echo 'OK'; // no test case for this function yet
?>
--EXPECT--
OK