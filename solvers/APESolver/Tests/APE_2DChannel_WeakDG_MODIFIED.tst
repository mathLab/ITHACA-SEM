<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D channel with base flow, P=5</description>
    <executable>APESolver</executable>
    <parameters>APE_2DChannel_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_2DChannel_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-12">0.102091</value>
            <value variable="u" tolerance="1e-12">0.000248662</value>
            <value variable="v" tolerance="1e-12">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-12">0.911223</value>
            <value variable="u" tolerance="1e-12">0.00221946</value>
            <value variable="v" tolerance="1e-12">3.54604e-17</value>
        </metric>
    </metrics>
</test>
