<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process npart to tecplot file </description>
    <executable python="true">chan3D_npart_tec.py</executable>
    <parameters></parameters>
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

