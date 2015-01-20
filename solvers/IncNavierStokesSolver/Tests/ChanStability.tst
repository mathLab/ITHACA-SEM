<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Linear stability (Mod. Arnoldi): Channel</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanStability.xml</parameters>
    <files>
        <file description="Session File">ChanStability.xml</file>
        <file description="Session File">ChanStability.bse</file>
        <file description="Session File">ChanStability.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.0372818</value>
            <value variable="v" tolerance="1e-12">0.0228428</value>
            <value variable="p" tolerance="1e-12">0.0181546</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.0272836</value>
            <value variable="v" tolerance="1e-12">0.0122607</value>
            <value variable="p" tolerance="1e-12">0.000130304</value>
        </metric>
    </metrics>
</test>
