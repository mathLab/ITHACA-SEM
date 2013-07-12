<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Courtemanche Cell model AF variant</description>
    <executable>PrePacing</executable>
    <parameters>CourtemancheAF.xml</parameters>
    <files>
        <file description="Session File">CourtemancheAF.xml</file>
    </files>
    <metrics>
        <metric type="Regex" id="1">
            <regex>
                ^#\s([\w]*)\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)
            </regex>
            <matches>
                <match>
                    <field>u</field>
                    <field tolerance="1e-06">-29.3514</field>
                </match>
                <match>
                    <field>m</field>
                    <field tolerance="1e-06">0.856243</field>
                </match>
                <match>
                    <field>h</field>
                    <field tolerance="1e-06">6.68791e-119</field>
                </match>
                <match>
                    <field>j</field>
                    <field tolerance="1e-06">1.61494e-10</field>
                </match>
                <match>
                    <field>o_a</field>
                    <field tolerance="1e-06">0.388397</field>
                </match>
                <match>
                    <field>o_i</field>
                    <field tolerance="1e-06">0.0394704</field>
                </match>
                <match>
                    <field>u_a</field>
                    <field tolerance="1e-06">0.548199</field>
                </match>
                <match>
                    <field>u_i</field>
                    <field tolerance="1e-06">0.994737</field>
                </match>
                <match>
                    <field>x_r</field>
                    <field tolerance="1e-06">0.119368</field>
                </match>
                <match>
                    <field>x_s</field>
                    <field tolerance="1e-06">0.070781</field>
                </match>
                <match>
                    <field>d</field>
                    <field tolerance="1e-06">0.0847185</field>
                </match>
                <match>
                    <field>f</field>
                    <field tolerance="1e-06">0.756132</field>
                </match>
                <match>
                    <field>f_Ca</field>
                    <field tolerance="1e-06">0.386251</field>
                </match>
                <match>
                    <field>U</field>
                    <field tolerance="1e-06">0.999994</field>
                </match>
                <match>
                    <field>V</field>
                    <field tolerance="1e-06">2.58224e-11</field>
                </match>
                <match>
                    <field>W</field>
                    <field tolerance="1e-06">0.983249</field>
                </match>
                <match>
                    <field>Na_i</field>
                    <field tolerance="1e-06">11.173</field>
                </match>
                <match>
                    <field>Ca_i</field>
                    <field tolerance="1e-06">0.000550987</field>
                </match>
                <match>
                    <field>K_i</field>
                    <field tolerance="1e-06">138.997</field>
                </match>
                <match>
                    <field>Ca_rel</field>
                    <field tolerance="1e-06">0.222938</field>
                </match>
                <match>
                    <field>Ca_up</field>
                    <field tolerance="1e-06">1.57489</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>




