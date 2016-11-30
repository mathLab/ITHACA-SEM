<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_m8.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_m8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">4.08607e-16</value>
            <value variable="v" tolerance="1e-12">1.54276e-16</value>
	    <value variable="p" tolerance="1e-12">1.1993e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">6.43929e-15</value>
            <value variable="v" tolerance="1e-12">6.36813e-16</value>
	    <value variable="p" tolerance="1e-12">4.86278e-14</value>
        </metric>
    </metrics>
</test>


