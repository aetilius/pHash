--TEST--
ph_audio_dist() function
--SKIPIF--
<?php 

if(!extension_loaded('pHash')) die('skip ');

if(!function_exists('ph_audio_dist')) die('skip not compiled in (HAVE_AUDIO_HASH)');

 ?>
--FILE--
<?php
echo 'OK'; // no test case for this function yet
?>
--EXPECT--
OK