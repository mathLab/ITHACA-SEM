<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>t106 step variant</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list 2d_bl_t106.mcf 2d_bl_t106-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">2d_bl_t106.mcf</file>
        <file description="Input File 2">2d_bl_t106.stp</file>
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
