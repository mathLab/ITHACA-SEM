<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Prismatic elements, P=6</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Prism_channel_m6_dealiasing.xml</parameters>
    <files>
        <file description="Session File">Prism_channel_m6_dealiasing.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">6.66372e-15</value>
            <value variable="v" tolerance="1e-12">7.15702e-15</value>
            <value variable="w" tolerance="1e-12">2.30265e-14</value>
	    <value variable="p" tolerance="1e-12">3.04193e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">3.1256e-14</value>
            <value variable="v" tolerance="1e-12">3.90951e-14</value>
            <value variable="w" tolerance="1e-12">7.31359e-14</value>
	    <value variable="p" tolerance="2e-12">1.93734e-12</value>
        </metric>
    </metrics>
</test>
