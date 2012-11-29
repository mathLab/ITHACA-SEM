<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=7 Body Force</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Channel_Flow_8modes_Body_Force.xml</parameters>
    <files>
        <file description="Session File">Channel_Flow_8modes_Body_Force.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">9.90432e-06</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-12">5.171e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.40068e-05</value>
            <value variable="v" tolerance="1e-12">4.82146e-16</value>
            <value variable="p" tolerance="1e-12">8.575e-14</value>
        </metric>
    </metrics>
</test>


