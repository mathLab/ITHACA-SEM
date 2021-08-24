<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=8 BodyForce</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_Standard_BodyForce.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_Standard_BodyForce.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">9.90432e-06</value>
            <value variable="v" tolerance="1e-12">2.55888e-16</value>
	    <value variable="p" tolerance="1e-12">1.25877e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.40068e-05</value>
            <value variable="v" tolerance="1e-12">8.04664e-16</value>
	    <value variable="p" tolerance="1e-12">7.56173e-14</value>
        </metric>
    </metrics>
</test>


