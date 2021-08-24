<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>.geo reader test mesh</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list 2d_bl.mcf 2d_bl-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">2d_bl.mcf</file>
        <file description="Input File 2">2d_bl.geo</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>.*Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
