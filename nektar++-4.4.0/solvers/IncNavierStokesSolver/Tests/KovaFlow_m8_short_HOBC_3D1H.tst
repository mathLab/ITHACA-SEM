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
            <value variable="u" tolerance="1e-11">2.53106e-13</value>
            <value variable="v" tolerance="1e-11">2.83952e-14</value>
            <value variable="w" tolerance="1e-11">7.00669e-14</value>
	    	<value variable="p" tolerance="1e-11">3.58285e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-11">2.56861e-12</value>
            <value variable="v" tolerance="1e-11">3.47405e-13</value>
            <value variable="w" tolerance="1e-11">5.99992e-13</value>
	    	<value variable="p" tolerance="1e-11">6.1374e-12</value>
        </metric>
    </metrics>
</test>
