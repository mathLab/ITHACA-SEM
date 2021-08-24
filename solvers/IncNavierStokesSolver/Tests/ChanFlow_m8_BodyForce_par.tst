<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=8 BodyForce</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>--use-scotch ChanFlow_m8_BodyForce.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">ChanFlow_m8_BodyForce.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">9.90432e-06</value>
            <value variable="v" tolerance="1e-8">2.82461e-11</value>
            <value variable="p" tolerance="1e-8">3.91926e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">1.40069e-05</value>
            <value variable="v" tolerance="1e-8">6.84991e-11</value>
            <value variable="p" tolerance="1e-8">1.21984e-08</value>
        </metric>
    </metrics>
</test>
