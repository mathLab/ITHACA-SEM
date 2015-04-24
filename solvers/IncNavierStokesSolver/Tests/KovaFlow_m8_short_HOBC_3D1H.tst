<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m8_short_HOBC_3D1H.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_m8_short_HOBC_3D1H.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-11">8.6117e-06</value>
            <value variable="v" tolerance="1e-11">3.04872e-14</value>
            <value variable="w" tolerance="1e-11">8.55267e-05</value>
	    	<value variable="p" tolerance="1e-11">9.07916e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-11">0.000369264</value>
            <value variable="v" tolerance="1e-11">4.00843e-13</value>
            <value variable="w" tolerance="1e-11">0.0022091</value>
	    	<value variable="p" tolerance="1e-11">0.00142775</value>
        </metric>
    </metrics>
</test>
