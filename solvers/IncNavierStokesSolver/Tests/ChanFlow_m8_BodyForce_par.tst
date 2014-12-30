<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=8 BodyForce</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_m8_BodyForce.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">ChanFlow_m8_BodyForce.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">9.9072e-06</value>
            <value variable="v" tolerance="1e-8">4.72839e-10</value>
            <value variable="p" tolerance="1e-8">1.19166e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">1.41125e-05</value>
            <value variable="v" tolerance="1e-8">1.5477e-09</value>
            <value variable="p" tolerance="1e-8">4.96317e-07</value>
        </metric>
    </metrics>
</test>
