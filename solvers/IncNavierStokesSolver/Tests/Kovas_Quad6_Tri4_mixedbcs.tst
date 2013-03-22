<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Steady Oseen Kovasznay flow, mixed elements and bcs P=7</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Kovas_Quad6_Tri4_mixedbcs.xml</parameters>
    <files>
        <file description="Session File">Kovas_Quad6_Tri4_mixedbcs.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.000529867</value>
            <value variable="v" tolerance="1e-12">0.000172531</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.00074911</value>
            <value variable="v" tolerance="1e-12">0.000353995</value>
        </metric>
    </metrics>
</test>
