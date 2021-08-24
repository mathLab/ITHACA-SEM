<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process npart to tecplot file </description>
    <executable python="true">Tet_channel_npart_wss.py</executable>
    <parameters></parameters>
    <files>
        <file description="Session File Directory">Tet_channel_m3_xml</file>
	<file description="Field File">Tet_channel_m3.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
           <value variable="x" tolerance="1e-6">0.5</value>
           <value variable="y" tolerance="1e-6">0</value>
           <value variable="z" tolerance="1e-6">0.288675</value>
           <value variable="Shear_x" tolerance="1e-6">2.06477e-16</value>
        </metric>
    </metrics>
</test>

