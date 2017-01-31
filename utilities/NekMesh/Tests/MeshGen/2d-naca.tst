<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Simple 2D NACA mesh</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list 2d-naca.mcf 2d-naca.xml:xml:test</parameters>
    <files>
        <file description="Input File">2d-naca.mcf</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
