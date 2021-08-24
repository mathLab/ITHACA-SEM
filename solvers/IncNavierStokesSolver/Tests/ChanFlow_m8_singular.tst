<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=8 Singularity Check</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_m8_singular.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_m8_singular.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">5.1664e-15</value>
            <value variable="v" tolerance="1e-6">5.99044e-15</value>
	    <value variable="p" tolerance="1e-6">4.74177e-12</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">2.74919e-14</value>
            <value variable="v" tolerance="1e-6">2.6359e-14</value>
	    <value variable="p" tolerance="1e-6">5.57478e-12</value>
        </metric>
    </metrics>
</test>


