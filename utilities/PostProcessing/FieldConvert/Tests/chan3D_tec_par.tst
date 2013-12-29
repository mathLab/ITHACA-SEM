<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3D tecplot output, par(2)</description>
    <executable>FieldConvert</executable>
    <parameters> chan3D.xml chan3D.fld chan3D.dat</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">chan3D.xml</file>
	<file description="Session File">chan3D.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="chan3D_P0.dat">
                <sha1>f9c494bd6582e356baa0feae344a79539ba4ce8d</sha1>
             </file>
         </metric>
        <metric type="file" id="2">
            <file filename="chan3D_P1.dat">
                <sha1>8e17ec2d41d7168848548a8f78e088b55cf1be05</sha1>
             </file>
         </metric>
    </metrics>
</test>

