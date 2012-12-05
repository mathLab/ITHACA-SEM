<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D channel flow, Tet elements, par(3), P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>-I GlobalSysSoln=IterativeMultiLevelStaticCond Tet_channel_m8_par.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Tet_channel_m8_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.98811e-14</value>
            <value variable="v" tolerance="1e-12">2.10102e-14</value>
            <value variable="w" tolerance="1e-12">1.09623e-13</value>
            <value variable="p" tolerance="1e-08">1.64394e-12</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">8.3757e-14</value>
            <value variable="v" tolerance="1e-12">9.56528e-14</value>
            <value variable="w" tolerance="1e-12">7.71383e-13</value>
            <value variable="p" tolerance="1e-08">7.46558e-12</value>
        </metric>
    </metrics>
</test>
