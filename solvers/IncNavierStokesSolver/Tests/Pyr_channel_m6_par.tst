<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D channel flow, Pyramidic elements, P=6</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>--use-scotch Pyr_channel_m6.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Pyr_channel_m6.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">3.77879e-14</value>
            <value variable="v" tolerance="1e-12">4.29463e-14</value>
            <value variable="w" tolerance="1e-11">3.96346e-13</value>
            <value variable="p" tolerance="1e-8">5.15486e-12</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">2.38835e-13</value>
            <value variable="v" tolerance="1e-12">2.2169e-13</value>
            <value variable="w" tolerance="1e-11">1.92338e-12</value>
            <value variable="p" tolerance="1e-8">2.83573e-11</value>
        </metric>
    </metrics>
</test>
