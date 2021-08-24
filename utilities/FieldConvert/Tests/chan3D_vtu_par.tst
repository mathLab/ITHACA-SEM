<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3D vtu output, par(2) </description>
    <executable>FieldConvert</executable>
    <parameters> -f chan3D.xml chan3D.fld chan3D.vtu</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">chan3D.xml</file>
	<file description="Session File">chan3D.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="chan3D_P0.vtu">
                <sha1>554d57eb11b057d88a03f0755f9922d7d3dfcb36</sha1>
             </file>
         </metric>
        <metric type="file" id="2">
            <file filename="chan3D_P1.vtu">
                <sha1>2b202e5e77c0cf184d62e4583bbb73fd069c0ca2</sha1>
             </file>
         </metric>
    </metrics>
</test>

