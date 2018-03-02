<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D channel flow, Hexahedral elements, P=3, flowrate</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>-I GlobalSysSoln=XxtStaticCond Hex_channel_m3_Flowrate.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Hex_channel_m3_Flowrate.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">5.94828e-16</value>
            <value variable="w" tolerance="1e-12">1.34946e-15</value>
            <value variable="p" tolerance="1e-8">9.43679e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">3.76927e-16</value>
            <value variable="v" tolerance="1e-12">1.70832e-15</value>
            <value variable="w" tolerance="1e-12">1.93179e-14</value>
            <value variable="p" tolerance="1e-8">1.49283e-13</value>
        </metric>
    </metrics>
</test>
