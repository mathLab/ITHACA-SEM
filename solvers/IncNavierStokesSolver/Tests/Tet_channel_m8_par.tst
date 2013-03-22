<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Tet elements, par(2), P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_channel_m8_par.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">Tet_channel_m8_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">1.04057e-11</value>
            <value variable="v" tolerance="1e-08">9.66528e-12</value>
            <value variable="w" tolerance="1e-08">1.09503e-10</value>
            <value variable="p" tolerance="1e-08">3.55555e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">6.08492e-11</value>
            <value variable="v" tolerance="1e-08">7.06681e-11</value>
            <value variable="w" tolerance="1e-08">2.66938e-09</value>
            <value variable="p" tolerance="1e-07">3.8533e-08</value>
        </metric>
    </metrics>
</test>
