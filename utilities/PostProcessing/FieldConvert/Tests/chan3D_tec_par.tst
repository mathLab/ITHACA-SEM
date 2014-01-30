<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3D tecplot output, par(2)</description>
    <executable>FieldConvert</executable>
    <parameters> -e chan3D.xml chan3D.fld chan3D.dat</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">chan3D.xml</file>
	<file description="Session File">chan3D.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="chan3D_P0.dat">
                <sha1>9bfd6d3bc24137f77f32a3148ee69eafdcf811ab</sha1>
             </file>
         </metric>
        <metric type="file" id="2">
            <file filename="chan3D_P1.dat">
                <sha1>7eb90cc0db9d011be2552d72a6c3c9e665a9d720</sha1>
             </file>
         </metric>
    </metrics>
</test>

