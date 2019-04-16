<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process npart to tecplot file </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e --nparts 2 chan3D_xml:xml chan3D.fld chan3D.plt </parameters>
    <files>
        <file description="Session File Directory">chan3D_xml</file>
	<file description="Field File">chan3D.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
           <value variable="x" tolerance="1e-6">1.17589</value>
           <value variable="y" tolerance="1e-6">1.1967</value>
           <value variable="z" tolerance="1e-6">1.21716</value>
           <value variable="u" tolerance="1e-6">1.49296</value>
        </metric>
    </metrics>
</test>

